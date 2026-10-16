module SpeedyTransformsReactantExt

using Reactant
using Reactant: Ops, TracedUtils
import AbstractFFTs
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays
import SpeedyWeatherInternals.Architectures: ReactantDevice, CPU, on_architecture
import SpeedyWeatherInternals.Utils: @maybe_jit, _jit

import SpeedyTransforms: SpectralTransform, ScratchMemory, ColumnScratchMemory, plan_FFTs,
    dense_legendre_polynomials, _transform_grid!, _transform_spec!, _transform_nonchunked!,
    _fourier!, _legendre!

# =====================================================================================
# FFT + LEGENDRE SPECTRAL TRANSFORM WITH REACTANT
#
# XLA compiles the whole transform (or the whole time step), so instead of many small operations
# (FFTW/cuFFT plans per ring and layer batch, Legendre loops or kernels with data-dependent loop
# bounds) the transform is written as few large array operations, a bit similar to the CUDA/HIP 
# graphs::
#
#  - Fourier: one `rfft`/`brfft` (`stablehlo.fft`) per run of latitudes with rings of equal length,
#    northern and southern rings together: one FFT per latitude on octahedral grids, one FFT for all
#    rings of a full grid. Only contiguous slices of the grid points, no gathers. No plans, XLA plans
#    the FFTs when it compiles.
#  - Legendre: one `dot_general` per direction, batched over the parity (even or odd degree l-m,
#    which gives the north/south symmetry) and the order m. The coefficients of order m (degrees
#    l = m..lmax, packed lower-triangular) are shifted to start at the first row and zero-padded
#    into a dense (parity, degree, order, layer) array, the polynomials are precomputed in the same
#    layout with the Legendre shortcut applied as zeros, see `dense_legendre_polynomials`.
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
Legendre polynomials in the dense (parity, degree, order, latitude ring) layout, flattened, zero
outside the lower triangle and where the Legendre shortcut `mmax_truncation` skips an order."""
function dense_legendre_polynomials(
        architecture::ReactantDevice,
        legendre_polynomials::LowerTriangularArray,
        mmax_truncation::AbstractVector{<:Integer},
    )
    (; lmax, mmax) = legendre_polynomials.spectrum              # 1-based
    polynomials = on_architecture(CPU(), legendre_polynomials.data)    # (lm, j)
    nlat_half = size(polynomials, 2)
    P = zeros(eltype(polynomials), 2, dense_lhalf(lmax), mmax, nlat_half)
    for j in 1:nlat_half
        lm = 0
        for m in 1:mmax, l in m:lmax
            lm += 1
            m <= mmax_truncation[j] + 1 || continue             # Legendre shortcut
            p, i = dense_position(l, m)
            P[p, i, m, j] = polynomials[lm, j]
        end
    end
    return on_architecture(architecture, vec(P))
end

"""$(TYPEDSIGNATURES)
Index maps between the packed lower-triangular coefficients (lm) and the dense (parity, degree,
order) layout, both as linear indices: `lm_of_dense` (length 2lhalf*mmax, `nharmonics + 1` for the
zero padding) and `dense_of_lm` (length nharmonics). Julia arrays, constants in the traced program."""
function dense_indices(spectrum::AbstractSpectrum)
    (; lmax, mmax) = spectrum
    lhalf = dense_lhalf(lmax)
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    lm_of_dense = fill(nharmonics + 1, 2, lhalf, mmax)
    dense_of_lm = zeros(Int, nharmonics)
    dense = LinearIndices(lm_of_dense)
    lm = 0
    for m in 1:mmax, l in m:lmax
        lm += 1
        p, i = dense_position(l, m)
        lm_of_dense[p, i, m] = lm
        dense_of_lm[lm] = dense[p, i, m]
    end
    return vec(lm_of_dense), dense_of_lm
end

# Fields and coefficients as (points or harmonics, columns) matrices: single-layer data are vectors
# (indexing them as matrices falls back to scalar indexing), views are materialized
as_matrix(data) = reshape(TracedUtils.materialize_traced_array(data), size(data, 1), :)
ncolumns(data) = length(data) ÷ size(data, 1)

# write the (n, columns) matrix `result` into `data` of any shape, or accumulate onto it
function write_matrix!(data, result, add::Bool)
    result = reshape(result, size(data))
    add ? (data .+= result) : (data .= result)
    return nothing
end

# real `dot_general` with a complex right-hand side: its real and imaginary parts are stacked as a
# new last dimension, which is the last dimension of the result, too
function dot_general_complex(real_lhs, complex_rhs; contracting_dimensions, batching_dimensions)
    stacked = cat(real.(complex_rhs), imag.(complex_rhs); dims = ndims(complex_rhs) + 1)
    result = Ops.dot_general(
        TracedUtils.materialize_traced_array(real_lhs), TracedUtils.materialize_traced_array(stacked);
        contracting_dimensions, batching_dimensions,
    )
    n = ndims(result)
    return complex.(selectdim(result, n, 1), selectdim(result, n, 2))
end

# =====================================================================================
# TRANSFORM: compile on the fly outside of compilation, Fourier and Legendre transform together
# =====================================================================================

function _transform_grid!(                  # SPECTRAL TO GRID
        field::AbstractField, coeffs::LowerTriangularArray,
        scratch_memory::ScratchMemory, S::ReactantSpectralTransform, unscale_coslat::Bool,
    )
    @maybe_jit S.architecture _transform_nonchunked!(field, coeffs, scratch_memory, S, unscale_coslat)
    return field
end

function _transform_spec!(                  # GRID TO SPECTRAL
        coeffs::LowerTriangularArray, field::AbstractField,
        scratch_memory::ScratchMemory, S::ReactantSpectralTransform,
    )
    @maybe_jit S.architecture _transform_nonchunked!(coeffs, field, scratch_memory, S)
    return coeffs
end

# =====================================================================================
# FOURIER
# =====================================================================================

"""$(TYPEDSIGNATURES)
Runs of latitudes `j1:j2` (1 to nlat_half, from the north pole to the equator) whose rings have the
same number of longitudes `nlon`, as pairs `nlon => j1:j2`. A ring and its southern counterpart
(ring nlat - j + 1) have the same length, so a run is Fourier-transformed in one FFT: one for a full
grid, one per latitude on octahedral grids. The northern rings of a run are contiguous grid points,
and so are its southern rings, which avoids gathers (whose reverse rules in Enzyme-MLIR fail for
the strided index patterns of north/south pairs)."""
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
(Forward) Fourier transform of all rings of `field` into the scratch memories `f_north`, `f_south`
(frequencies beyond a ring's `nlon ÷ 2 + 1` zero, as is the southern equator ring), see `_fourier!`.
One FFT per run of latitudes with rings of equal length, see `latitude_runs`."""
function _fourier!(                         # GRID TO SPECTRAL
        f_north::AbstractArray{<:Any, 3},
        f_south::AbstractArray{<:Any, 3},
        field::AbstractField,
        S::ReactantSpectralTransform,
    )
    (; rings, nlons, nlat, nfreq_max) = S
    (; nlat_half) = S.grid
    data = as_matrix(field.data)
    K = size(data, 2)

    NF = Reactant.unwrapped_eltype(eltype(field))      # Float32 for a traced field of TracedRNumber{Float32}
    @assert NF == eltype(S) "Number format of grid $NF and SpectralTransform $(eltype(S)) need to match."
    @boundscheck K <= size(f_north, 2) || throw(DimensionMismatch(S, field))

    north_spectra, south_spectra = [], []   # per run: (nfreq_max, latitudes, K)
    for (nlon, js) in latitude_runs(nlons)
        n = length(js)
        # northern rings j1:j2, then the southern ones (latitudes reversed), (nlon, 2n, K)
        points = cat(data[northern_points(rings, js), :], data[southern_points(rings, js, nlat), :]; dims = 1)
        spectrum = AbstractFFTs.rfft(reshape(points, nlon, 2n, K), 1)                     # (nfreq, 2n, K)
        nfreq = size(spectrum, 1)
        if nfreq < nfreq_max                # zero-pad to the longest ring
            spectrum = cat(spectrum, zeros(eltype(spectrum), nfreq_max - nfreq, 2n, K); dims = 1)
        end
        push!(north_spectra, spectrum[:, 1:n, :])
        push!(south_spectra, n > 1 ? reverse(spectrum[:, (n + 1):2n, :]; dims = 2) : spectrum[:, (n + 1):2n, :])
    end
    north = cat(north_spectra...; dims = 2) # (nfreq_max, nlat_half, K)
    south = cat(south_spectra...; dims = 2)

    if isodd(nlat)                          # the equator ring is only transformed once, as north
        equator_mask = reshape([j < nlat_half for j in 1:nlat_half], 1, nlat_half, 1)
        south = south .* equator_mask
    end

    view(f_north, :, 1:K, :) .= permutedims(north, (1, 3, 2))
    view(f_south, :, 1:K, :) .= permutedims(south, (1, 3, 2))
    return nothing
end

"""$(TYPEDSIGNATURES)
Inverse Fourier transform of the Legendre-transformed `g_north`, `g_south` into `field` (or
accumulated with `add = true`), see `_fourier!`. One FFT per run of latitudes with rings of equal
length, see `latitude_runs`."""
function _fourier!(                         # SPECTRAL TO GRID
        field::AbstractField,
        g_north::AbstractArray{<:Any, 3},
        g_south::AbstractArray{<:Any, 3},
        S::ReactantSpectralTransform;
        add::Bool = false,
    )
    (; nlons, nlat) = S
    K = ncolumns(field.data)

    @boundscheck K <= size(g_north, 2) || throw(DimensionMismatch(S, field))

    north = permutedims(g_north[:, 1:K, :], (1, 3, 2))      # (nfreq_max, nlat_half, K)
    south = permutedims(g_south[:, 1:K, :], (1, 3, 2))

    north_points, south_points = [], []     # per run: (nlon * latitudes, K) in grid order
    runs = latitude_runs(nlons)
    for (nlon, js) in runs
        nfreq = nlon ÷ 2 + 1
        n = length(js)
        south_spectra = n > 1 ? reverse(south[1:nfreq, js, :]; dims = 2) : south[1:nfreq, js, :]   # grid order
        spectra = cat(north[1:nfreq, js, :], south_spectra; dims = 2)                   # (nfreq, 2n, K)
        points = reshape(AbstractFFTs.brfft(spectra, nlon, 1), nlon * 2n, K)            # northern, southern rings
        push!(north_points, points[1:(nlon * n), :])
        push!(south_points, points[(nlon * n + 1):end, :])
    end
    if isodd(nlat)                          # drop the southern copy of the equator ring (last run, first ring)
        nlon_equator = first(runs[end])
        south_points[end] = south_points[end][(nlon_equator + 1):end, :]
    end
    # northern hemisphere from the pole, then southern hemisphere from the equator (runs in reverse)
    all_points = cat(north_points..., reverse(south_points)...; dims = 1)
    write_matrix!(field.data, all_points, add)
    return nothing
end

# =====================================================================================
# LEGENDRE
# =====================================================================================

"""$(TYPEDSIGNATURES)
Inverse Legendre transform of `specs` into `g_north`, `g_south` (frequencies beyond `mmax` zero)
as one batched matrix multiply over parity and order, see `_legendre!`."""
function _legendre!(                        # SPECTRAL TO GRID
        g_north::AbstractArray{<:Any, 3},
        g_south::AbstractArray{<:Any, 3},
        specs::LowerTriangularArray,
        scratch_memory::ColumnScratchMemory,  # unused, for the CPU transform only
        S::ReactantSpectralTransform;
        unscale_coslat::Bool = false,
    )
    (; spectrum, lon_offsets, coslat⁻¹, nfreq_max) = S
    (; lmax, mmax) = spectrum
    (; nlat_half) = S.grid
    lhalf = dense_lhalf(lmax)
    data = as_matrix(specs.data)
    K = size(data, 2)

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck K <= size(g_north, 2) || throw(DimensionMismatch(S, specs))

    # coefficients into the dense (parity, degree, order, layer) layout, zero-padded
    lm_of_dense, _ = dense_indices(spectrum)
    padded = cat(data, zeros(eltype(data), 1, K); dims = 1)
    coeffs = reshape(padded[lm_of_dense, :], 2, lhalf, mmax, K)
    polynomials = reshape(S.legendre_polynomials_dense, 2, lhalf, mmax, nlat_half)

    # P[p, i, m, j] * coeffs[p, i, m, k] summed over the degree i -> [p, m, j, k]
    sums = dot_general_complex(polynomials, coeffs; batching_dimensions = ([1, 3], [1, 3]), contracting_dimensions = ([2], [2]))
    even = sums[1, :, :, :]                 # (mmax, nlat_half, K)
    odd = sums[2, :, :, :]

    # longitude offset rotation (and unscaling by coslat) per order and ring
    rotation = unscale_coslat ? lon_offsets .* reshape(coslat⁻¹[1:nlat_half], 1, nlat_half) : lon_offsets
    north = rotation .* (even .+ odd)       # north = even + odd, south = even - odd
    south = rotation .* (even .- odd)

    zero_frequencies = zeros(eltype(north), nfreq_max - mmax, nlat_half, K)
    view(g_north, :, 1:K, :) .= permutedims(cat(north, zero_frequencies; dims = 1), (1, 3, 2))
    view(g_south, :, 1:K, :) .= permutedims(cat(south, zero_frequencies; dims = 1), (1, 3, 2))
    return nothing
end

"""$(TYPEDSIGNATURES)
(Forward) Legendre transform of the Fourier-transformed `f_north`, `f_south` into `specs` (or
accumulated with `add = true`) as one batched matrix multiply over parity and order, see `_legendre!`."""
function _legendre!(                        # GRID TO SPECTRAL
        specs::LowerTriangularArray,
        f_north::AbstractArray{<:Any, 3},
        f_south::AbstractArray{<:Any, 3},
        scratch_memory::ColumnScratchMemory,  # unused, for the CPU transform only
        S::ReactantSpectralTransform;
        add::Bool = false,
    )
    (; spectrum, lon_offsets, solid_angles) = S
    (; lmax, mmax) = spectrum
    (; nlat_half) = S.grid
    lhalf = dense_lhalf(lmax)
    K = ncolumns(specs.data)

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck K <= size(f_north, 2) || throw(DimensionMismatch(S, specs))

    fn = f_north[1:mmax, 1:K, :]            # (mmax, K, nlat_half)
    fs = f_south[1:mmax, 1:K, :]

    # solid angle quadrature weights and longitude offset, rotated back to the prime meridian
    weights = reshape(solid_angles[1:nlat_half], 1, nlat_half) .* conj.(lon_offsets)
    weights = reshape(weights, mmax, 1, nlat_half)
    even = reshape(weights .* (fn .+ fs), 1, mmax, K, nlat_half)
    odd = reshape(weights .* (fn .- fs), 1, mmax, K, nlat_half)
    symmetric = cat(even, odd; dims = 1)    # (parity, order, layer, ring)
    polynomials = reshape(S.legendre_polynomials_dense, 2, lhalf, mmax, nlat_half)

    # P[p, i, m, j] * symmetric[p, m, k, j] summed over the rings j -> [p, m, i, k]
    sums = dot_general_complex(polynomials, symmetric; batching_dimensions = ([1, 3], [1, 2]), contracting_dimensions = ([4], [4]))
    dense = reshape(permutedims(sums, (1, 3, 2, 4)), 2 * lhalf * mmax, K)

    _, dense_of_lm = dense_indices(spectrum)
    coeffs = dense[dense_of_lm, :]
    write_matrix!(specs.data, coeffs, add)
    return nothing
end

end
