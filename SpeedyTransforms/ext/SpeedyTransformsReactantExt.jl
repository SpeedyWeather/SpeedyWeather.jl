module SpeedyTransformsReactantExt

using Reactant
using Reactant: Ops, TracedUtils
import AbstractFFTs
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays
import SpeedyWeatherInternals.Architectures: ReactantDevice, CPU, on_architecture
import SpeedyWeatherInternals.Utils: @maybe_jit

import SpeedyTransforms: SpectralTransform, ScratchMemory, plan_FFTs, dense_legendre_polynomials,
    dense_fourier_matrices, _transform_grid!, _transform_spec!

# =====================================================================================
# SPECTRAL TRANSFORM WITH REACTANT
#
# XLA compiles the whole transform (or the whole time step), so instead of many small operations
# (FFTW/cuFFT plans per ring and layer batch, Legendre loops or kernels with data-dependent loop
# bounds) the transform is written as few large array operations, a bit similar to the CUDA/HIP
# graphs. Fourier and Legendre transform are one traced function that passes their intermediate as a
# value (XLA allocates it and chooses its layout) instead of through scratch memory, and real throughout
# (real and imaginary parts as a dimension of their own) until the packed spectral coefficients:
#
#   grid → spectral: Fourier → (ring, order, real/imaginary, north/south, K) → Legendre → coefficients
#   spectral → grid: coefficients → Legendre → (order, ring, real/imaginary, north/south, K) → Fourier
#
# (the two layouts differ as `dot_general` puts the batch dimension of its result first).
#
#  - Fourier, two ways (`fourier_dft` of the `SpectralTransform`), with the same interface: forward the
#    orders m = 0..mmax-1 of every ring times the ring's solid angle (quadrature weight) and rotated back
#    by its longitude offset, inverse rotated by the longitude offset first
#    - as truncated DFT (`fourier_dft = true`, the default up to `DFT_MAX_TRUNCATION`): one batched
#      matrix multiply per direction over all rings (zero-padded to the longest), see "FOURIER AS
#      TRUNCATED DFT";
#    - with FFTs (`fourier_dft = false`): one `rfft`/`brfft` (`stablehlo.fft`) per run of latitudes
#      with rings of equal length, northern and southern rings together: one FFT per latitude on
#      octahedral grids, one FFT for all rings of a full grid. No plans, XLA plans the FFTs when it
#      compiles. XLA can't record FFTs into its CUDA graphs though, so on GPUs this splits a compiled
#      program into several graphs. See "FOURIER WITH FFTs".
#  - Legendre: one `dot_general` per direction, batched over the parity (even or odd degree l-m,
#    which gives the north/south symmetry) and the order m. The coefficients of order m (degrees
#    l = m..lmax, packed lower-triangular) are shifted to start at the first row and zero-padded
#    into a dense (parity, order, degree) array, the layout of the forward transform's `dot_general`
#    result; the polynomials are precomputed in the same layout with the Legendre shortcut applied
#    as zeros, see `dense_indices` and `dense_legendre_polynomials`.
#
# Everything here is traced: outside of compilation (e.g. for initial conditions) the transform is
# compiled on the fly with `@maybe_jit`. The CUDA graphs of SpeedyTransformsCUDAExt are never used,
# XLA runs a compiled program as one CUDA graph itself.
# =====================================================================================

const ReactantSpectralTransform = SpectralTransform{NF, <:ReactantDevice} where {NF}

# No FFT plans on Reactant, XLA plans the FFTs (`stablehlo.fft`) when compiling
plan_FFTs(
    ::ReactantDevice, planned_K::AbstractVector{<:Integer}, ::Type, grid::AbstractGrid, nlayers::Integer,
    scratch_memory_north::AbstractArray, rings, nlons::Vector{<:Int},
) = (nothing, nothing, nothing, nothing)

# DENSE LAYOUT, see the comment above. Degree index l' = l - m + 1 within the column of order m
# (1-based), split by parity p = 1 (l' odd, 0-based even degree l-m) or 2 (l' even) into
# l' = 2(i - 1) + p, i = 1..lhalf with lhalf = ceil(lmax / 2).
dense_position(l::Integer, m::Integer) = (lprime = l - m + 1; (2 - isodd(lprime), (lprime + 1) ÷ 2))
dense_lhalf(lmax::Integer) = (lmax + 1) ÷ 2

"""$(TYPEDSIGNATURES)
Index maps between the packed lower-triangular coefficients (lm) and the dense (parity, order,
degree) layout, both as linear indices: `lm_of_dense` (length 2mmax*lhalf, `nharmonics + 1` for the
zero padding) and `dense_of_lm` (length nharmonics). Julia arrays, constants in the traced program."""
function dense_indices(spectrum::AbstractSpectrum)
    (; lmax, mmax) = spectrum
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    lm_of_dense = fill(nharmonics + 1, 2, mmax, dense_lhalf(lmax))
    dense_of_lm = zeros(Int, nharmonics)
    dense = LinearIndices(lm_of_dense)
    lm = 0
    for m in 1:mmax, l in m:lmax
        lm += 1
        p, i = dense_position(l, m)
        lm_of_dense[p, m, i] = lm
        dense_of_lm[lm] = dense[p, m, i]
    end
    return vec(lm_of_dense), dense_of_lm
end

"""$(TYPEDSIGNATURES)
Legendre polynomials in the dense (parity, order, degree, latitude ring) layout, flattened, zero
outside the lower triangle and where the Legendre shortcut `mmax_truncation` skips an order."""
function dense_legendre_polynomials(
        architecture::ReactantDevice,
        legendre_polynomials::LowerTriangularArray,
        mmax_truncation::AbstractVector{<:Integer},
    )
    (; spectrum) = legendre_polynomials
    (; lmax, mmax) = spectrum                                   # 1-based
    polynomials = on_architecture(CPU(), legendre_polynomials.data)    # (lm, j)
    nlat_half = size(polynomials, 2)
    _, dense_of_lm = dense_indices(spectrum)
    P = zeros(eltype(polynomials), 2, mmax, dense_lhalf(lmax), nlat_half)
    reshape(P, :, nlat_half)[dense_of_lm, :] = polynomials
    for j in 1:nlat_half                    # Legendre shortcut: orders m > mmax_truncation[j] (0-based) of ring j
        P[:, (mmax_truncation[j] + 2):mmax, :, j] .= 0
    end
    return on_architecture(architecture, vec(P))
end

materialize(x) = TracedUtils.materialize_traced_array(x)

# `a` and `b` of the same size stacked along a new dimension `dims` (of size 2), e.g. the real and
# imaginary parts, the northern and southern rings
stack_pair(a, b, dims::Integer) = cat(insert_dimension(a, dims), insert_dimension(b, dims); dims)
insert_dimension(a, dims::Integer) = reshape(a, size(a)[1:(dims - 1)]..., 1, size(a)[dims:end]...)

# Fields and coefficients as (points or harmonics, columns) matrices: single-layer data are vectors
# (indexing them as matrices falls back to scalar indexing), views are materialized
as_matrix(data) = reshape(materialize(data), size(data, 1), :)

# write the (n, columns) matrix `result` into `data` of any shape
write_matrix!(data, result) = (data .= reshape(result, size(data)); nothing)

# =====================================================================================
# TRANSFORM: compile on the fly outside of compilation, Fourier and Legendre transform together
# =====================================================================================

function _transform_grid!(                  # SPECTRAL TO GRID
        field::AbstractField, coeffs::LowerTriangularArray,
        scratch_memory::ScratchMemory,      # unused, the transform needs no scratch memory on Reactant
        S::ReactantSpectralTransform, unscale_coslat::Bool,
    )
    @maybe_jit S.architecture _transform_reactant!(field, coeffs, S, unscale_coslat)
    return field
end

function _transform_spec!(                  # GRID TO SPECTRAL
        coeffs::LowerTriangularArray, field::AbstractField,
        scratch_memory::ScratchMemory,      # unused, the transform needs no scratch memory on Reactant
        S::ReactantSpectralTransform,
    )
    @maybe_jit S.architecture _transform_reactant!(coeffs, field, S)
    return coeffs
end

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid) on Reactant: inverse Legendre transform, then the inverse
Fourier transform as truncated DFT (`S.fourier_dft`) or with FFTs, see the comment at the top."""
function _transform_reactant!(              # SPECTRAL TO GRID
        field::AbstractField, coeffs::LowerTriangularArray, S::ReactantSpectralTransform, unscale_coslat::Bool,
    )
    @boundscheck SpeedyTransforms.ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck SpeedyTransforms.ismatching(S, coeffs) || throw(DimensionMismatch(S, coeffs))
    NF = Reactant.unwrapped_eltype(eltype(field))      # Float32 for a traced field of TracedRNumber{Float32}
    @assert NF == eltype(S) "Number format of grid $NF and SpectralTransform $(eltype(S)) need to match."

    spectra = legendre_inverse(as_matrix(coeffs.data), S, unscale_coslat)
    data = S.fourier_dft ? dft_fourier_inverse(spectra, S) : fft_fourier_inverse(spectra, S)
    write_matrix!(field.data, data)
    return field
end

"""$(TYPEDSIGNATURES)
Spectral transform (grid to spectral) on Reactant: Fourier transform as truncated DFT
(`S.fourier_dft`) or with FFTs, then the Legendre transform, see the comment at the top."""
function _transform_reactant!(              # GRID TO SPECTRAL
        coeffs::LowerTriangularArray, field::AbstractField, S::ReactantSpectralTransform,
    )
    @boundscheck SpeedyTransforms.ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck SpeedyTransforms.ismatching(S, coeffs) || throw(DimensionMismatch(S, coeffs))
    NF = Reactant.unwrapped_eltype(eltype(field))
    @assert NF == eltype(S) "Number format of grid $NF and SpectralTransform $(eltype(S)) need to match."

    data = as_matrix(field.data)
    spectra = S.fourier_dft ? dft_fourier_forward(data, S) : fft_fourier_forward(data, S)
    write_matrix!(coeffs.data, legendre_forward(spectra, S))
    return coeffs
end

# =====================================================================================
# LEGENDRE
# =====================================================================================

"""$(TYPEDSIGNATURES)
Forward Legendre transform from the weighted and rotated Fourier coefficients `spectra` (nlat_half,
mmax, real/imaginary, north/south, K) of `dft_fourier_forward` or `fft_fourier_forward` to the packed
spectral coefficients (nharmonics, K): one matrix multiply batched over parity and order with the real
and imaginary parts as a free dimension."""
function legendre_forward(spectra, S::ReactantSpectralTransform)
    (; spectrum) = S
    (; lmax, mmax) = spectrum
    (; nlat_half) = S.grid
    lhalf = dense_lhalf(lmax)
    K = size(spectra, 5)

    north = spectra[:, :, :, 1, :]          # (nlat_half, mmax, real/imaginary, K)
    south = spectra[:, :, :, 2, :]
    symmetric = stack_pair(north .+ south, north .- south, 1)  # (parity: even/odd, ring, order, real/imaginary, K)
    polynomials = reshape(S.legendre_polynomials_dense, 2, mmax, lhalf, nlat_half)

    # P[p, m, i, j] * symmetric[p, j, m, c, k] summed over the rings j -> [p, m, i, c, k]
    sums = Ops.dot_general(
        materialize(polynomials), materialize(symmetric);
        batching_dimensions = ([1, 2], [1, 3]), contracting_dimensions = ([4], [2]),
    )
    _, dense_of_lm = dense_indices(spectrum)
    parts = reshape(sums, 2 * mmax * lhalf, 2, K)[dense_of_lm, :, :]
    return complex.(parts[:, 1, :], parts[:, 2, :])
end

"""$(TYPEDSIGNATURES)
Inverse Legendre transform from the packed spectral coefficients `data` (nharmonics, K) to the real
and imaginary parts of the Fourier coefficients of all rings (mmax, nlat_half, real/imaginary,
north/south, K), unscaled by the cosine of latitude with `unscale_coslat`, not yet rotated by the
longitude offsets (that's part of the inverse Fourier transform): one matrix multiply batched over
parity and order with the real and imaginary parts as a free dimension."""
function legendre_inverse(data, S::ReactantSpectralTransform, unscale_coslat::Bool)
    (; spectrum, coslat⁻¹) = S
    (; lmax, mmax) = spectrum
    (; nlat_half) = S.grid
    lhalf = dense_lhalf(lmax)
    K = size(data, 2)

    # coefficients into the dense (parity, order, degree, real/imaginary, K) layout, zero-padded
    lm_of_dense, _ = dense_indices(spectrum)
    parts = stack_pair(real.(data), imag.(data), 2)                 # (nharmonics, real/imaginary, K)
    padded = cat(parts, zeros(eltype(parts), 1, 2, K); dims = 1)
    coeffs = reshape(padded[lm_of_dense, :, :], 2, mmax, lhalf, 2, K)
    polynomials = reshape(S.legendre_polynomials_dense, 2, mmax, lhalf, nlat_half)

    # P[p, m, i, j] * coeffs[p, m, i, c, k] summed over the degree i -> [p, m, j, c, k]
    sums = Ops.dot_general(
        materialize(polynomials), materialize(coeffs);
        batching_dimensions = ([1, 2], [1, 2]), contracting_dimensions = ([3], [3]),
    )
    even = sums[1, :, :, :, :]              # (mmax, nlat_half, real/imaginary, K)
    odd = sums[2, :, :, :, :]
    north = even .+ odd                     # north = even + odd, south = even - odd
    south = even .- odd
    if unscale_coslat
        coslat⁻¹_north = reshape(coslat⁻¹[1:nlat_half], 1, nlat_half, 1, 1)
        north = north .* coslat⁻¹_north
        south = south .* coslat⁻¹_north
    end
    return stack_pair(north, south, 4)      # (mmax, nlat_half, real/imaginary, north/south, K)
end

# =====================================================================================
# FOURIER WITH FFTs (fourier_dft = false)
# =====================================================================================

"""$(TYPEDSIGNATURES)
Runs of latitudes `j1:j2` (1 to nlat_half, from the north pole to the equator) whose rings have the
same number of longitudes `nlon`, as pairs `nlon => j1:j2`. A ring and its southern counterpart
(ring nlat - j + 1) have the same length, so a run is Fourier-transformed in one FFT: one for a full
grid, one per latitude on octahedral grids. The northern rings of a run are contiguous grid points,
and so are its southern rings, which avoids gathers (whose reverse rules in Enzyme-MLIR failed for
the strided index patterns of north/south pairs before Reactant 0.2.294)."""
function latitude_runs(nlons::AbstractVector{<:Integer})
    runs = Pair{Int, UnitRange{Int}}[]
    for (j, nlon) in enumerate(nlons)
        if !isempty(runs) && first(runs[end]) == nlon
            runs[end] = nlon => first(last(runs[end])):j
        else
            push!(runs, nlon => j:j)
        end
    end
    return runs
end

# grid points of the northern rings j1:j2 and of their southern counterparts (rings nlat-j2+1:nlat-j1+1)
northern_points(rings, js) = first(rings[first(js)]):last(rings[last(js)])
southern_points(rings, js, nlat) = first(rings[nlat - last(js) + 1]):last(rings[nlat - first(js) + 1])

"""$(TYPEDSIGNATURES)
Fourier transform of all rings of `data` (npoints, K) with one `rfft` per run of latitudes with rings
of equal length (see `latitude_runs`), as `dft_fourier_forward`: the orders m = 0..mmax-1 (zero beyond a
ring's Nyquist frequency, and for the southern copy of the equator ring with an odd number of rings),
times the solid angle of the ring and rotated back by its longitude offset,
(nlat_half, mmax, real/imaginary, north/south, K)."""
function fft_fourier_forward(data, S::ReactantSpectralTransform)
    (; rings, nlons, nlat, lon_offsets, solid_angles) = S
    (; nlat_half) = S.grid
    (; mmax) = S.spectrum
    K = size(data, 2)

    north_spectra, south_spectra = [], []   # per run: (mmax, latitudes, K)
    for (nlon, js) in latitude_runs(nlons)
        n = length(js)
        # northern rings j1:j2, then the southern ones (latitudes reversed), (nlon, 2n, K)
        points = cat(data[northern_points(rings, js), :], data[southern_points(rings, js, nlat), :]; dims = 1)
        spectrum = AbstractFFTs.rfft(reshape(points, nlon, 2n, K), 1)                     # (nfreq, 2n, K)
        nfreq = size(spectrum, 1)
        spectrum = nfreq >= mmax ? spectrum[1:mmax, :, :] :                              # orders 0..mmax-1
            cat(spectrum, zeros(eltype(spectrum), mmax - nfreq, 2n, K); dims = 1)        # zero beyond the ring
        push!(north_spectra, spectrum[:, 1:n, :])
        push!(south_spectra, n > 1 ? reverse(spectrum[:, (n + 1):2n, :]; dims = 2) : spectrum[:, (n + 1):2n, :])
    end
    north = cat(north_spectra...; dims = 2) # (mmax, nlat_half, K)
    south = cat(south_spectra...; dims = 2)

    if isodd(nlat)                          # the equator ring is only transformed once, as north
        equator_mask = reshape([j < nlat_half for j in 1:nlat_half], 1, nlat_half, 1)
        south = south .* equator_mask
    end

    # solid angle quadrature weights and longitude offset, rotated back to the prime meridian,
    # into the (ring, order, real/imaginary) layout of the Legendre transform
    weights = reshape(solid_angles[1:nlat_half], 1, nlat_half) .* conj.(lon_offsets)   # (mmax, nlat_half)
    function real_parts(spectra)
        weighted = permutedims(weights .* spectra, (2, 1, 3))                              # (nlat_half, mmax, K)
        return stack_pair(real.(weighted), imag.(weighted), 3)
    end
    return stack_pair(real_parts(north), real_parts(south), 4)
end

"""$(TYPEDSIGNATURES)
Inverse Fourier transform of the real and imaginary parts `spectra` (mmax, nlat_half, real/imaginary,
north/south, K) of the orders m = 0..mmax-1 of all rings to the grid points (npoints, K), as
`dft_fourier_inverse`: rotated by the longitude offset, then one `brfft` per run of latitudes with
rings of equal length (see `latitude_runs`)."""
function fft_fourier_inverse(spectra, S::ReactantSpectralTransform)
    (; nlons, nlat, nfreq_max, lon_offsets) = S
    (; mmax) = S.spectrum
    K = size(spectra, 5)

    # complex coefficients rotated by the longitude offset, (mmax, nlat_half, K), zero-padded to the
    # frequencies of the longest ring
    north = lon_offsets .* complex.(spectra[:, :, 1, 1, :], spectra[:, :, 2, 1, :])
    south = lon_offsets .* complex.(spectra[:, :, 1, 2, :], spectra[:, :, 2, 2, :])
    if nfreq_max > mmax
        zero_frequencies = zeros(eltype(north), nfreq_max - mmax, size(north, 2), K)
        north = cat(north, zero_frequencies; dims = 1)
        south = cat(south, zero_frequencies; dims = 1)
    end

    north_points, south_points = [], []     # per run: (nlon * latitudes, K) in grid order
    runs = latitude_runs(nlons)
    for (nlon, js) in runs
        nfreq = nlon ÷ 2 + 1
        n = length(js)
        south_spectra = n > 1 ? reverse(south[1:nfreq, js, :]; dims = 2) : south[1:nfreq, js, :]   # grid order
        spectra_run = cat(north[1:nfreq, js, :], south_spectra; dims = 2)               # (nfreq, 2n, K)
        points = reshape(AbstractFFTs.brfft(spectra_run, nlon, 1), nlon * 2n, K)        # northern, southern rings
        push!(north_points, points[1:(nlon * n), :])
        push!(south_points, points[(nlon * n + 1):end, :])
    end
    if isodd(nlat)                          # drop the southern copy of the equator ring (last run, first ring)
        nlon_equator = first(runs[end])
        if size(south_points[end], 1) == nlon_equator
            pop!(south_points)              # a run of the equator ring only (e.g. octahedral grids)
        else
            south_points[end] = south_points[end][(nlon_equator + 1):end, :]
        end
    end
    # northern hemisphere from the pole, then southern hemisphere from the equator (runs in reverse)
    return cat(north_points..., reverse(south_points)...; dims = 1)
end

# =====================================================================================
# FOURIER AS TRUNCATED DFT (fourier_dft = true)
#
# Only the first mmax of a ring's nlon/2 + 1 Fourier frequencies are used, so the Fourier transform of a
# ring is a truncated DFT, a small dense matrix. With all rings zero-padded to the longest one, the
# Fourier transform of all rings is one matrix multiply per direction, batched over the latitude rings,
# which XLA keeps in its CUDA graph (FFTs it can't record into one). The quadrature weights and the
# longitude offsets are folded into the matrices:
#
#   grid → spectral: gather the rings (nlon_max, nlat_half, north/south, K) → DFT → (ring, order,
#       real/imaginary, north/south, K)
#   spectral → grid: (order, ring, real/imaginary, north/south, K) → inverse DFT → (ring, longitude,
#       north/south, K) → gather the grid points
#
# The index tables of the gathers are Julia arrays, constants in the traced program.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Forward and inverse matrices of the Fourier transform as a truncated DFT, in Float64 on the host, for
rings with `nlons` longitudes (northern hemisphere, a southern ring has the matrices of its northern
counterpart), their longitude offsets `lon_offsets` (mmax, nlat_half), solid angles `solid_angles` and
the orders m = 0..mmax-1.
- `forward` (2mmax, nlon_max, nlat_half) maps the grid points of a ring (zero-padded to `nlon_max`) to the
  real (rows 1:mmax) and imaginary parts (rows mmax+1:2mmax) of its Fourier coefficients as `rfft` does,
  times the ring's solid angle and rotated back by its longitude offset. Orders beyond the ring's Nyquist
  frequency are zero.
- `inverse` (nlon_max, 2mmax, nlat_half) maps the real (columns 1:mmax) and imaginary parts (columns
  mmax+1:2mmax) of the coefficients, first rotated by the longitude offset, to the grid points of a ring
  as `brfft` does: factor 2 for 0 < m < nlon/2 (the other half of the spectrum), the imaginary part of
  m = 0 and of the Nyquist frequency (after the rotation) ignored, zero beyond the ring.
The phases `exp(±2πi m (i-1) / nlon)` are computed with the argument reduced exactly in integers."""
function fourier_matrices(
        nlons::AbstractVector{<:Integer},
        lon_offsets::AbstractMatrix,
        solid_angles::AbstractVector,
        mmax::Integer,
    )
    nlat_half = length(nlons)
    nlon_max = maximum(nlons)
    forward = zeros(Float64, 2mmax, nlon_max, nlat_half)
    inverse = zeros(Float64, nlon_max, 2mmax, nlat_half)
    for (j, nlon) in enumerate(nlons)
        for m in 0:min(mmax - 1, nlon ÷ 2)          # orders up to the ring's Nyquist frequency
            offset = ComplexF64(lon_offsets[m + 1, j])
            weight = Float64(solid_angles[j]) * conj(offset)
            factor = (m == 0 || 2m == nlon) ? 1 : 2 # brfft: the other half of the spectrum by symmetry
            for i in 1:nlon
                phase = cispi(2 * mod(m * (i - 1), nlon) / nlon)    # exp(2πi m (i-1) / nlon)
                f = weight * conj(phase)
                forward[m + 1, i, j] = real(f)
                forward[mmax + m + 1, i, j] = imag(f)
                b = factor * offset * phase         # Re(c * b) = Re(c) Re(b) - Im(c) Im(b)
                inverse[i, m + 1, j] = real(b)
                inverse[i, mmax + m + 1, j] = -imag(b)
            end
        end
    end
    return forward, inverse
end

"""$(TYPEDSIGNATURES)
The matrices of `fourier_matrices` in number format `NF` on Reactant, flattened, if `fourier_dft`,
otherwise empty vectors."""
function dense_fourier_matrices(
        architecture::ReactantDevice, ::Type{NF}, fourier_dft::Bool, nlons, lon_offsets, solid_angles, mmax,
    ) where {NF}
    fourier_dft || return on_architecture(architecture, zeros(NF, 0)), on_architecture(architecture, zeros(NF, 0))
    forward, inverse = fourier_matrices(nlons, lon_offsets, solid_angles, mmax)
    return on_architecture(architecture, vec(NF.(forward))), on_architecture(architecture, vec(NF.(inverse)))
end

"""$(TYPEDSIGNATURES)
Index tables of the DFT as linear indices, Julia arrays (constants in the traced program):
- `gather` (length nlon_max * nlat_half * 2): grid point at longitude i of the northern (1) or southern
  (2) ring j in a (nlon_max, nlat_half, 2) array, or npoints + 1 (a zero row appended to the grid points)
  for the zero padding beyond the ring and for the southern copy of the equator ring of an odd number of
  rings (transformed once, as northern ring).
- `positions` (length npoints): position of every grid point in the (nlat_half, nlon_max, 2) result of
  the inverse DFT, the equator ring from its northern copy."""
function ring_indices(S::ReactantSpectralTransform)
    (; rings, nlat, nlon_max) = S
    (; nlat_half) = S.grid
    npoints = sum(length, rings)
    gather = fill(npoints + 1, nlon_max, nlat_half, 2)
    positions = zeros(Int, npoints)
    for j in 1:nlat_half, hemisphere in 1:2
        ring = hemisphere == 1 ? j : nlat - j + 1
        hemisphere == 2 && ring == j && continue    # equator ring (odd nlat): northern copy only
        for (i, point) in enumerate(rings[ring])
            gather[i, j, hemisphere] = point
            positions[point] = j + (i - 1) * nlat_half + (hemisphere - 1) * nlat_half * nlon_max
        end
    end
    return vec(gather), positions
end

"""$(TYPEDSIGNATURES)
Fourier transform of all rings of `data` (npoints, K) as one matrix multiply with the truncated DFT
matrix (see `fourier_matrices`), batched over the rings: the orders m = 0..mmax-1, times the solid angle
of the ring and rotated back by its longitude offset, (nlat_half, mmax, real/imaginary, north/south, K)."""
function dft_fourier_forward(data, S::ReactantSpectralTransform)
    (; nlon_max) = S
    (; nlat_half) = S.grid
    (; mmax) = S.spectrum
    K = size(data, 2)

    gather, _ = ring_indices(S)
    padded = cat(data, zeros(eltype(data), 1, K); dims = 1)                 # zero row for the padding
    rings = reshape(padded[gather, :], nlon_max, nlat_half, 2, K)
    matrix = reshape(S.fourier_matrix_forward, 2mmax, nlon_max, nlat_half)

    # W[r, i, j] * rings[i, j, h, k] summed over the longitudes i -> [j, r, h, k]
    spectra = Ops.dot_general(
        materialize(matrix), materialize(rings);
        batching_dimensions = ([3], [2]), contracting_dimensions = ([2], [1]),
    )
    return reshape(spectra, nlat_half, mmax, 2, 2, K)
end

"""$(TYPEDSIGNATURES)
Inverse Fourier transform of the real and imaginary parts `spectra` (mmax, nlat_half, real/imaginary,
north/south, K) of the orders m = 0..mmax-1 of all rings to the grid points (npoints, K), one matrix
multiply with the truncated inverse DFT matrix (see `fourier_matrices`), batched over the rings."""
function dft_fourier_inverse(spectra, S::ReactantSpectralTransform)
    (; nlon_max) = S
    (; nlat_half) = S.grid
    (; mmax) = S.spectrum
    K = size(spectra, 5)

    _, positions = ring_indices(S)
    matrix = reshape(S.fourier_matrix_inverse, nlon_max, mmax, 2, nlat_half)

    # B[i, m, c, j] * spectra[m, j, c, h, k] summed over orders m and real/imaginary parts c -> [j, i, h, k]
    rings = Ops.dot_general(
        materialize(matrix), materialize(spectra);
        batching_dimensions = ([4], [2]), contracting_dimensions = ([2, 3], [1, 3]),
    )
    return reshape(rings, nlat_half * nlon_max * 2, K)[positions, :]
end

end
