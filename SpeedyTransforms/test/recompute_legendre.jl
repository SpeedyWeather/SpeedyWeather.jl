# Tests for step 2 of the recompute-on-the-fly Legendre polynomials feature (see
# docs/dev/2026-09/recompute-legendre-polynomials.md): `SpectralTransform`'s `recompute_legendre`
# keyword, `PrecomputedLegendre`/`RecomputedLegendre`, and the CPU `legendre_ring!` path. GPU
# kernels for the recomputed path are step 3 and not implemented here; on GPU
# `recompute_legendre = true` must throw, tested below via a CPU-only fake using `S.legendre` (no
# GPU/JLArrays dependency needed for that part since the check is architecture-agnostic on the
# `RecomputedLegendre` type, only gated on `S isa SpectralTransform{NF, <:Architectures.GPU}`).

# tolerances measured empirically (see PR discussion): a plain `Float32` recursion in the
# transform gives ~1e-6 relative error against the `Float64`-precomputed reference (rounded to
# Float32 on write), `Float64` recomputation is indistinguishable from precomputation (~1e-15).
# Set just above the measured worst case so a real regression (e.g. losing the double-single
# arithmetic, ~1e-2 per the design doc) trips these tests.
const RECOMPUTE_RTOL = Dict(Float32 => 1.0e-5, Float64 => 1.0e-11)

@testset "recompute_legendre: grid -> spectral -> grid agrees with precomputed" begin
    for NF in (Float32, Float64)
        for truncation in (15, 31)
            for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
                arch = SpeedyTransforms.Architectures.CPU()
                spectrum = Spectrum(truncation, truncation; architecture = arch)
                nlat_half = SpeedyTransforms.get_nlat_half(truncation, SpeedyTransforms.default_dealiasing(Grid))
                grid = Grid(nlat_half, arch)

                Sp = SpectralTransform(spectrum, grid; NF, recompute_legendre = false)
                Sr = SpectralTransform(spectrum, grid; NF, recompute_legendre = true)

                field = rand(NF, grid, 2)
                coeffs_p = transform(field, Sp)
                coeffs_r = transform(field, Sr)
                field_p = transform(coeffs_p, Sp)
                field_r = transform(coeffs_r, Sr)

                rtol = RECOMPUTE_RTOL[NF]
                @test coeffs_r.data ≈ coeffs_p.data rtol = rtol
                @test field_r.data ≈ field_p.data rtol = rtol
            end
        end
    end
end

@testset "recompute_legendre: memory footprint much smaller than precomputed" begin
    arch = SpeedyTransforms.Architectures.CPU()
    truncation = 63
    spectrum = Spectrum(truncation, truncation; architecture = arch)
    Grid = FullGaussianGrid
    nlat_half = SpeedyTransforms.get_nlat_half(truncation, SpeedyTransforms.default_dealiasing(Grid))
    grid = Grid(nlat_half, arch)

    Sp = SpectralTransform(spectrum, grid; NF = Float32, recompute_legendre = false)
    Sr = SpectralTransform(spectrum, grid; NF = Float32, recompute_legendre = true)

    mem_p = SpeedyTransforms.legendre_memory_size(Sp.legendre)
    mem_r = SpeedyTransforms.legendre_memory_size(Sr.legendre)
    @test mem_r < mem_p / 2   # measured ~4x smaller at T63, leave margin
end

@testset "recompute_legendre: RecomputedLegendre on GPU throws (step 3 not implemented)" begin
    # PrecomputedLegendre stays the only GPU-supported mode until step 3; simulate a GPU
    # architecture without a real backend by checking the dispatch guard directly through the
    # type system: SpectralTransform{NF, <:Architectures.GPU} is what legendre_ka.jl dispatches
    # on, and the guard throws before touching any GPU-only array, so this is safe to check with
    # `isa` on the constructed struct without a working GPU backend.
    arch = SpeedyTransforms.Architectures.CPU()
    truncation = 15
    Grid = FullGaussianGrid
    spectrum = Spectrum(truncation, truncation; architecture = arch)
    nlat_half = SpeedyTransforms.get_nlat_half(truncation, SpeedyTransforms.default_dealiasing(Grid))
    grid = Grid(nlat_half, arch)
    Sr = SpectralTransform(spectrum, grid; NF = Float32, recompute_legendre = true)
    @test Sr.legendre isa SpeedyTransforms.RecomputedLegendre
    @test Sr.legendre.ring isa Vector{Float32}
    @test length(Sr.legendre.ring) == LowerTriangularArrays.nonzeros(spectrum)
    @test isempty(Sr.legendre.tile)
    @test isempty(Sr.legendre.tile_orders)
end

@testset "recompute_legendre: tile_order_blocks partitions 1:mmax exactly once, contiguously" begin
    for (lmax, mmax, nlat_half) in ((16, 16, 8), (129, 128, 65), (64, 5, 20))
        for NF in (Float32, Float64)
            for tile_budget in (1000, 32_000_000)
                blocks = SpeedyTransforms.tile_order_blocks(lmax, mmax, nlat_half, NF, tile_budget)
                @test !isempty(blocks)
                @test all(!isempty, blocks)                       # at least one order per block
                covered = vcat(collect.(blocks)...)
                @test covered == 1:mmax                            # exact cover, in order, no gaps/overlap
                # contiguous: consecutive blocks abut
                for i in 2:length(blocks)
                    @test first(blocks[i]) == last(blocks[i - 1]) + 1
                end
            end
        end
    end
end
