"""Storage for the associated Legendre polynomials used by `SpectralTransform`, either
precomputed once and stored (`PrecomputedLegendre`, today's behaviour) or recomputed on the fly
from the `ScaledLegendre` recursion (`RecomputedLegendre`), trading `O(nonzeros(spectrum) *
nlat_half)` of memory for `O(nonzeros(spectrum) + nlat_half * mmax)` at the cost of running the
recursion inside the transform. See `docs/dev/2026-09/recompute-legendre-polynomials.md` for the
full design. `SpectralTransform` holds one `AbstractLegendrePolynomials` in its `legendre` field,
chosen at construction time by the `recompute_legendre` keyword; every reader downstream goes
through `legendre_ring!` (CPU, `legendre.jl`) or reads the `PrecomputedLegendre` fields directly
(GPU, `legendre_ka.jl` — the GPU recomputed kernels are a later step)."""
abstract type AbstractLegendrePolynomials{NF} end

"""$(TYPEDSIGNATURES)
Legendre polynomials precomputed once for every spherical harmonic `(l, m)` and latitude ring `j`
(northern hemisphere only, by symmetry), exactly as `SpectralTransform` has always done. Fields
$(TYPEDFIELDS)"""
struct PrecomputedLegendre{NF, LowerTriangularArrayType, VectorType} <: AbstractLegendrePolynomials{NF}
    "Legendre polynomials λ_l^m(cos(colat_j)), stored as (lm, j)"
    polynomials::LowerTriangularArrayType
    "Transposed copy (j, lm), flattened so it needs no type parameter of its own; index it as
    [(lm - 1) * nlat_half + j]. Only the GPU inverse transform reads this layout (the forward
    transform reads `polynomials` directly); empty on CPU."
    polynomials_transposed::VectorType
end

Adapt.@adapt_structure PrecomputedLegendre

function Architectures.on_architecture(arch::AbstractArchitecture, L::PrecomputedLegendre{NF}) where {NF}
    polynomials = on_architecture(arch, L.polynomials)
    polynomials_transposed = on_architecture(arch, L.polynomials_transposed)
    return PrecomputedLegendre{NF, typeof(polynomials), typeof(polynomials_transposed)}(
        polynomials, polynomials_transposed,
    )
end

"""$(TYPEDSIGNATURES)
Precompute the Legendre polynomials for `spectrum`/`grid` in number format `NF`, moved onto
`architecture`. This is exactly the code `SpectralTransform`'s constructor used to run inline
(moved here verbatim so the precomputed path stays bit-for-bit identical): loop over latitude
rings `j`, run `AssociatedLegendrePolynomials.jl`'s `Legendre.λlm!` once per ring, and — on GPU
only — keep a second, transposed copy for the inverse transform kernel (see the struct docs)."""
function PrecomputedLegendre(
        ::Type{NF}, spectrum::AbstractSpectrum, grid::AbstractGrid,
        cos_colat::AbstractVector, architecture::AbstractArchitecture,
    ) where {NF}
    (; lmax, mmax) = spectrum                             # 1-based max degree l, order m
    (; nlat_half) = grid

    polynomials = zeros(LowerTriangularArray{NF}, spectrum, nlat_half)
    polynomials_j = zeros(NF, lmax, mmax)                 # temporary for one latitude
    for j in 1:nlat_half                                  # only one hemisphere due to symmetry
        Legendre.λlm!(polynomials_j, lmax - 1, mmax - 1, cos_colat[j])   # precompute l, m 0-based
        polynomials[:, j] = LowerTriangularArray(polynomials_j)         # store
    end

    # TRANSPOSED COPY OF THE LEGENDRE POLYNOMIALS FOR THE GPU INVERSE TRANSFORM, see the struct
    # field. Costs as much memory again as the polynomials themselves, so only built on GPU.
    polynomials_transposed = architecture isa Architectures.GPU ?
        vec(permutedims(polynomials.data)) :
        on_architecture(architecture, zeros(NF, 0))

    return PrecomputedLegendre{NF, typeof(polynomials), typeof(polynomials_transposed)}(
        polynomials, polynomials_transposed,
    )
end

"""$(TYPEDSIGNATURES)
Legendre polynomials recomputed on the fly from the `ScaledLegendre` recursion instead of stored.
Holds only `O(nonzeros(spectrum) + nlat_half * mmax)` of coefficients/starting values rather than
the `O(nonzeros(spectrum) * nlat_half)` full table. Fields $(TYPEDFIELDS)

`ring` and `tile` are mutually-exclusive scratch buffers for the two architectures: the CPU path
(`legendre_ring!` in `legendre.jl`) recomputes one whole ring (all columns, all orders) at a time
into `ring`, while the GPU forward-transform path (step 3) recomputes one *tile* — a contiguous
block of orders `m`, see `tile_orders` — at a time into `tile`. Only one of the two is ever
non-empty for a given architecture."""
struct RecomputedLegendre{NF, VectorType, MatrixType, IntMatrixType, RingType} <: AbstractLegendrePolynomials{NF}
    "1-based max degree l of the spectrum this was built for; needed to recover each order's
    column length (lmax - m + 1) since it is not otherwise stored on this struct."
    lmax::Int
    "Recursion coefficients α, β (double-single hi/lo pairs), length nonzeros(spectrum) + 2,
    indexed by the running lm index, see `ScaledLegendre.recursion_coefficients`."
    αhi::VectorType
    αlo::VectorType
    βhi::VectorType
    βlo::VectorType
    "Sectoral (diagonal, l == m) starting values λ_m^m(cos_colat[j]) and their extended-exponent
    scale, (nlat_half, mmax), see `ScaledLegendre.sectoral_modes`."
    sectoral_hi::MatrixType
    sectoral_lo::MatrixType
    sectoral_scale::IntMatrixType
    "cos(colat) as a double-single pair, length nlat_half each, in number format NF. Split rather
    than stored as a single NF value because it multiplies the running state at *every* recursion
    step: rounding it to Float32 alone would inject a relative error of eps(Float32) that the
    double-single arithmetic downstream cannot recover. `xlo` is exactly zero for NF == Float64."
    xhi::VectorType
    xlo::VectorType
    "CPU-only scratch: one recomputed ring (all columns) at a time, length nonzeros(spectrum) on
    CPU, length 0 on GPU."
    ring::RingType
    "GPU-only scratch for the forward transform (step 3): one tile (block of orders) at a time,
    (nnz_tile, nlat_half); length 0 on CPU."
    tile::MatrixType
    "Blocks of orders m that partition 1:mmax for the tiled GPU forward transform, host side
    (plain Vector, never moved to the device); empty on CPU. Each block is a contiguous range of
    orders and therefore also a contiguous range of the running lm index (columns are stored
    order-by-order in the lower-triangular layout), so every coefficient lm belongs to exactly one
    block and the tiled forward transform writes each coefficient exactly once."
    tile_orders::Vector{UnitRange{Int}}
end

Adapt.@adapt_structure RecomputedLegendre

function Architectures.on_architecture(arch::AbstractArchitecture, L::RecomputedLegendre{NF}) where {NF}
    αhi = on_architecture(arch, L.αhi)
    αlo = on_architecture(arch, L.αlo)
    βhi = on_architecture(arch, L.βhi)
    βlo = on_architecture(arch, L.βlo)
    sectoral_hi = on_architecture(arch, L.sectoral_hi)
    sectoral_lo = on_architecture(arch, L.sectoral_lo)
    sectoral_scale = on_architecture(arch, L.sectoral_scale)
    xhi = on_architecture(arch, L.xhi)
    xlo = on_architecture(arch, L.xlo)
    tile = on_architecture(arch, L.tile)
    return RecomputedLegendre{NF, typeof(αhi), typeof(sectoral_hi), typeof(sectoral_scale), typeof(L.ring)}(
        L.lmax, αhi, αlo, βhi, βlo, sectoral_hi, sectoral_lo, sectoral_scale,
        xhi, xlo, L.ring, tile, L.tile_orders,
    )
end

"""$(TYPEDSIGNATURES)
Partition the orders `1:mmax` of a spectrum with `lmax` degrees into consecutive blocks ("tiles")
for the GPU forward transform, such that each block's total column count (summed `lmax - m + 1`
over the orders `m` in the block) times `nlat_half * sizeof(NF)` bytes stays under `tile_budget`,
with at least one order per block (so a single very wide order is never split). Blocks are
contiguous ranges of `m`, and therefore also contiguous ranges of the running `lm` index (see
`RecomputedLegendre`'s `tile_orders` field docs) — that is what lets the tiled forward transform
write every coefficient exactly once, with no accumulation across tiles."""
function tile_order_blocks(lmax::Integer, mmax::Integer, nlat_half::Integer, ::Type{NF}, tile_budget::Integer) where {NF}
    max_cols = max(1, tile_budget ÷ (nlat_half * sizeof(NF)))   # max total lm-rows per tile
    blocks = UnitRange{Int}[]
    m_start = 1
    cols = 0
    for m in 1:mmax
        column_length = lmax - m + 1
        if cols > 0 && cols + column_length > max_cols
            push!(blocks, m_start:(m - 1))
            m_start = m
            cols = 0
        end
        cols += column_length
    end
    push!(blocks, m_start:mmax)
    return blocks
end

"""$(TYPEDSIGNATURES)
Build a `RecomputedLegendre` for `spectrum`/`grid` in number format `NF`, moved onto
`architecture`. The recursion coefficients and sectoral starting values are always built on the
CPU (in `Float64`, only split into the `NF` double-single representation at the end, see
`ScaledLegendre`) and then transferred with `on_architecture`. `tile_budget` (bytes, default 32
MB) sizes the GPU forward-transform scratch `tile`, see `tile_order_blocks`; unused on CPU."""
function RecomputedLegendre(
        ::Type{NF}, spectrum::AbstractSpectrum, grid::AbstractGrid,
        cos_colat::AbstractVector, architecture::AbstractArchitecture;
        tile_budget::Integer = 32_000_000,
    ) where {NF}
    (; lmax, mmax) = spectrum                              # 1-based max degree l, order m
    (; nlat_half) = grid
    n = LowerTriangularArrays.nonzeros(spectrum)

    # recursion coefficients and sectoral starting values, always computed in Float64 on the CPU
    # and only split into the NF double-single representation at the end, see ScaledLegendre.
    αhi, αlo, βhi, βlo = ScaledLegendre.recursion_coefficients(NF, lmax, mmax)
    sectoral_hi, sectoral_lo, sectoral_scale = ScaledLegendre.sectoral_modes(NF, cos_colat, mmax)
    xhi = zeros(NF, nlat_half)
    xlo = zeros(NF, nlat_half)
    for j in 1:nlat_half        # split through Float64, see the xhi/xlo field docs
        xhi[j], xlo[j] = ScaledLegendre.split_two_float(NF, Float64(cos_colat[j]))
    end

    on_gpu = architecture isa Architectures.GPU

    # `ring` (CPU) and `tile` (GPU) are mutually exclusive scratch buffers, see the struct docs.
    ring = zeros(NF, on_gpu ? 0 : n)
    if on_gpu
        tile_orders = tile_order_blocks(lmax, mmax, nlat_half, NF, tile_budget)
        max_block_columns = maximum(sum(lmax - m + 1 for m in block) for block in tile_orders)
        tile = zeros(NF, max_block_columns, nlat_half)
    else
        tile_orders = UnitRange{Int}[]
        tile = zeros(NF, 0, 0)
    end

    αhi_d = on_architecture(architecture, αhi)
    αlo_d = on_architecture(architecture, αlo)
    βhi_d = on_architecture(architecture, βhi)
    βlo_d = on_architecture(architecture, βlo)
    sectoral_hi_d = on_architecture(architecture, sectoral_hi)
    sectoral_lo_d = on_architecture(architecture, sectoral_lo)
    sectoral_scale_d = on_architecture(architecture, sectoral_scale)
    xhi_d = on_architecture(architecture, xhi)
    xlo_d = on_architecture(architecture, xlo)
    tile_d = on_architecture(architecture, tile)

    return RecomputedLegendre{
        NF, typeof(αhi_d), typeof(sectoral_hi_d), typeof(sectoral_scale_d), typeof(ring),
    }(
        lmax,
        αhi_d, αlo_d, βhi_d, βlo_d,
        sectoral_hi_d, sectoral_lo_d, sectoral_scale_d,
        xhi_d, xlo_d, ring, tile_d, tile_orders,
    )
end

"""$(TYPEDSIGNATURES)
`Base.summarysize`-based memory footprint of an `AbstractLegendrePolynomials`, in bytes, summing
every array field (the plain-`Int`/`Vector{UnitRange{Int}}` bookkeeping fields are negligible next
to the coefficient/polynomial arrays and are included via `summarysize` for simplicity)."""
legendre_memory_size(L::AbstractLegendrePolynomials) = Base.summarysize(L)

"""$(TYPEDSIGNATURES)
Short label for the mode an `AbstractLegendrePolynomials` was built in, used by `show`."""
legendre_mode(::PrecomputedLegendre) = "precomputed"
legendre_mode(::RecomputedLegendre) = "recomputed"
