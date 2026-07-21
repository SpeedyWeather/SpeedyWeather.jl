# Function barrier for batched or serial transforms. FFTW/cuFFT plans bake the batch dim K
# into the plan, so we look up a pre-planned bundle by K = size(field, 2) and fall back to
# the serial path (K=1 plan, looped) when no batched plan exists for that K. The batch dim 
# K is the sum of number of vertical layers batched together in a single FFT plan.
function _fourier!(f_north, f_south, field::AbstractField, S::SpectralTransform)
    K = size(field, 2)
    if K > 1 && haskey(S.rfft_plans_batched, K)
        return _fourier_batched!(f_north, f_south, field, S)
    else
        return _fourier_serial!(f_north, f_south, field, S)
    end
end

function _fourier!(field::AbstractField, f_north, f_south, S::SpectralTransform; add::Bool = false)
    K = size(field, 2)
    if K > 1 && haskey(S.brfft_plans_batched, K)
        return _fourier_batched!(field, f_north, f_south, S; add)
    else
        return _fourier_serial!(field, f_north, f_south, S; add)
    end
end

"""$(TYPEDSIGNATURES)
(Forward) FFT, applied in zonal direction of `field` provided. Depending on
the field's array type (CPU/GPU) the actual application of the FFT is slightly
different because of strided array/view support of the GPU FFT libraries."""
function _apply_batched_fft!(
        f_out::AbstractArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
        j::Int,
        nfreq::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true
    )
    nlayers = size(field, 2)        # number of vertical layers
    rfft_plan = S.rfft_plans_batched[nlayers][j]  # concrete K-batched plan for this ring

    # Perform the FFT
    if not_equator # skip FFT, redundant when north already did that latitude
        view(f_out, 1:nfreq, 1:nlayers, j) .= rfft_plan * view_only_on_cpu(field.data, ilons, :)
    else
        f_out[1:nfreq, 1:nlayers, j] .= 0
    end
    return nothing
end

# CPU version with view
@inline view_only_on_cpu(A::AbstractArray, I...) = view(A, I...)
# GPU version without view
@inline view_only_on_cpu(A::AbstractGPUArray, I...) = A[I...]

"""$(TYPEDSIGNATURES)
(Inverse) FFT, applied in zonal direction of `field` provided. Depending on
the field's array type (CPU/GPU) the actual application of the FFT is slightly
different because of strided array/view support of the GPU FFT libraries."""
function _apply_batched_fft!(
        field::AbstractField,
        g_in::AbstractArray{<:Complex, 3},
        S::SpectralTransform,
        j::Int,
        nlon::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true,
        add::Bool = false,          # accumulate onto `field` instead of overwriting? (used by Enzyme adjoint rule)
    )
    nlayers = size(field, 2)        # number of vertical layers
    brfft_plan = S.brfft_plans_batched[nlayers][j]  # concrete K-batched plan for this ring
    nfreq = nlon ÷ 2 + 1

    if not_equator  # skip FFT, redundant when north already did that latitude
        dest = view(field.data, ilons, :)
        rhs = brfft_plan * view_only_on_cpu(g_in, 1:nfreq, 1:nlayers, j)
        add ? (dest .+= rhs) : (dest .= rhs)
    end
    return nothing
end

# Forward FFT Application Function for serial FFT (K=1 plans, applied per layer)
function _apply_serial_fft!(
        f_out::AbstractArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
        j::Int,
        k::Int,
        nfreq::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true
    )
    rfft_plan = S.rfft_plan_serial[j]   # concrete K=1 (1-D) plan, planned wrt nlon on ring
    k_grid = eachlayer(field)[k]     # Precomputed ring index (as a Cartesian index)

    if not_equator
        view(f_out, 1:nfreq, k, j) .= rfft_plan * view(field.data, ilons, k_grid)
    else
        fill!(view(f_out, 1:nfreq, k, j), 0)
    end
    return nothing
end

# Inverse FFT Application Function for serial FFT (K=1 plans, applied per layer)
function _apply_serial_fft!(
        field::AbstractField,
        g_in::AbstractArray{<:Complex, 3},
        S::SpectralTransform,
        j::Int,
        k::Int,
        nfreq::Int,
        ilons::UnitRange{Int};
        not_equator::Bool = true,
        add::Bool = false,          # accumulate onto `field` instead of overwriting? (used by Enzyme adjoint rule)
    )
    brfft_plan = S.brfft_plan_serial[j]   # concrete K=1 (1-D) plan, planned wrt nlon on ring
    k_grid = eachlayer(field)[k]       # Precomputed ring index (as a Cartesian index)

    if not_equator
        dest = view(field.data, ilons, k_grid)
        rhs = brfft_plan * view(g_in, 1:nfreq, k, j)
        add ? (dest .+= rhs) : (dest .= rhs)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
(Forward) Fast Fourier transform (grid to spectral) in zonal direction of `field`,
stored in scratch memories `f_north`, `f_south` to be passed on to the Legendre transform.
Batched version that requires the number of vertical layers to be the same as precomputed in `S`.
Not to be called directly, use `transform!` instead."""
function _fourier_batched!(                 # GRID TO SPECTRAL
        f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
        f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
        field::AbstractField,                   # gridded input
        S::SpectralTransform                   # precomputed transform
    )
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck haskey(S.rfft_plans_batched, nlayers) || throw(DimensionMismatch(S, field))
    # scratch dim 2 is the per-call capacity (= max(planned_K) on CPU, nlayers elsewhere);
    # allow it to exceed nlayers so the bound passes for both full-K and chunked calls.
    @boundscheck (size(f_north) == size(f_south) && size(f_north, 1) == S.nfreq_max &&
                  size(f_north, 3) == nlat_half && size(f_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, field))

    return @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # corresponding southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring
        nfreq = nlon ÷ 2 + 1                  # linear max Fourier frequency wrt to nlon
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        ilons = rings[j_north]              # in-ring indices northern ring
        # FOURIER TRANSFORM in zonal direction, northern latitude
        _apply_batched_fft!(f_north, field, S, j, nfreq, ilons)

        # and southern latitude if not on Equator
        ilons = rings[j_south]              # in-ring indices southern ring
        _apply_batched_fft!(f_south, field, S, j, nfreq, ilons; not_equator = not_equator)

    end
end

"""$(TYPEDSIGNATURES)
(Forward) Fast Fourier transform (grid to spectral) in zonal direction of `field`,
stored in scratch memories `f_north`, `f_south` to be passed on to the Legendre transform.
Serial version that does not require the number of vertical layers to be the same as precomputed in `S`.
Not to be called directly, use `transform!` instead."""
function _fourier_serial!(                  # GRID TO SPECTRAL
        f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
        f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
        field::AbstractField,                   # gridded input
        S::SpectralTransform                   # precomputed transform
    )
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))
    @boundscheck (size(f_north) == size(f_south) && size(f_north, 1) == S.nfreq_max &&
                  size(f_north, 3) == nlat_half && size(f_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, field))

    return @inbounds for (k, k_grid) in zip(1:nlayers, eachlayer(field))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j = j_north                         # symmetric index / ring-away from pole index
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j]                     # number of longitudes on this ring (north or south)
            nfreq = nlon ÷ 2 + 1                 # linear max Fourier frequency wrt to nlon
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            ilons = rings[j_north]              # in-ring indices northern ring

            # Apply FFT in the northern latitudes
            _apply_serial_fft!(f_north, field, S, j, k, nfreq, ilons)

            # southern latitude, don't call redundant 2nd fft if ring is on equator
            ilons = rings[j_south]                      # in-ring indices southern ring
            _apply_serial_fft!(f_south, field, S, j, k, nfreq, ilons; not_equator = not_equator)
        end
    end
end

"""$(TYPEDSIGNATURES)
Inverse fast Fourier transform (spectral to grid) of Legendre-transformed inputs `g_north` and `g_south`
to be stored in `field`. Not to be called directly, use `transform!` instead."""
function _fourier_batched!(                 # SPECTRAL TO GRID
        field::AbstractField,                   # gridded output
        g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
        g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
        S::SpectralTransform;                  # precomputed transform
        add::Bool = false,                      # accumulate onto `field`? (used by Enzyme adjoint rule)
    )
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck haskey(S.brfft_plans_batched, nlayers) || throw(DimensionMismatch(S, field))   # otherwise FFTW complains
    @boundscheck (size(g_north) == size(g_south) && size(g_north, 1) == S.nfreq_max &&
                  size(g_north, 3) == nlat_half && size(g_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, field))

    return @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring (north or south)
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        ilons = rings[j_north]              # in-ring indices northern ring

        # northern latitude
        _apply_batched_fft!(field, g_north, S, j, nlon, ilons; add = add)

        # southern latitude, don't call redundant 2nd FFT if ring is on equator
        ilons = rings[j_south]              # in-ring indices southern ring
        _apply_batched_fft!(field, g_south, S, j, nlon, ilons; not_equator = not_equator, add = add)

    end
end

"""$(TYPEDSIGNATURES)
(Inverse) Fast Fourier transform (spectral to grid) of Legendre-transformed inputs `g_north` and `g_south`
to be stored in `field`. Serial version that does not require the number of vertical layers to be the same
as precomputed in `S`. Not to be called directly, use `transform!` instead."""
function _fourier_serial!(                  # SPECTRAL TO GRID
        field::AbstractField,                   # gridded output
        g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
        g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
        S::SpectralTransform;                  # precomputed transform
        add::Bool = false,                      # accumulate onto `field`? (used by Enzyme adjoint rule)
    )
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))     # otherwise FFTW complains
    @boundscheck (size(g_north) == size(g_south) && size(g_north, 1) == S.nfreq_max &&
                  size(g_north, 3) == nlat_half && size(g_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, field))

    return @inbounds for (k, k_grid) in zip(1:nlayers, eachlayer(field))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j = j_north                         # symmetric index / ring-away from pole index
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j]                     # number of longitudes on this ring (north or south)
            nfreq = nlon ÷ 2 + 1                 # linear max Fourier frequency wrt to nlon
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            # Apply FFT in the northern latitudes
            ilons = rings[j_north]              # in-ring indices northern ring
            _apply_serial_fft!(field, g_north, S, j, k, nfreq, ilons; add = add)

            # southern latitude, don't call redundant 2nd fft if ring is on equator
            ilons = rings[j_south]              # in-ring indices southern ring
            _apply_serial_fft!(field, g_south, S, j, k, nfreq, ilons; not_equator = not_equator, add = add)
        end
    end
end

"""$(TYPEDSIGNATURES)
Generate the per-ring FFT plans. FFTW/cuFFT plans bake the batch dim K into the plan at
construction, so we build **separate, concretely-typed** plan sets (they are stored on `S` with
concrete types so `plan * view` / `mul!` is a static call — see `SpectralTransform` fields):

- The **serial (K=1)** plans (1-D input, one layer) are always built — the per-layer fallback used by
  `_fourier_serial!`. Returned as concretely-typed `Vector`s.
- The **batched (K>1)** plans (2-D input, batched over the trailing dim) are built for each K in
  `planned_K` and returned as `Dict{Int, Vector{P}}` keyed by K, with a concrete plan type `P`. When no
  batched K is planned (e.g. `nlayers==1`) the Dict is empty (never indexed).

Because the K=1 (1-D) and K>1 (2-D) plans are different concrete FFTW/cuFFT types, they cannot share a
single `Dict`, hence the serial/batched split. `fake_grid_data` must have `size(_, 2) ≥ maximum(planned_K)`
so the per-K views below are valid (in practice sized to max K = `S.nlayers`)."""
function plan_FFTs(
        planned_K::AbstractVector{<:Integer},
        fake_grid_data::AbstractField{NF, N, <:AbstractArray{NF}},
        scratch_memory_north::AbstractArray{Complex{NF}},
        rings,
        nlons::Vector{<:Int}
    ) where {NF <: AbstractFloat, N}

    nlat_half = length(nlons)
    nfreqj(j) = nlons[j] ÷ 2 + 1

    # SERIAL (K=1, 1-D) plans — always built. Comprehensions give concretely-typed `Vector`s.
    rfft_plan_serial = [AbstractFFTs.plan_rfft(view_only_on_cpu(fake_grid_data.data, rings[j], 1), 1)
                        for j in 1:nlat_half]
    brfft_plan_serial = [AbstractFFTs.plan_brfft(view_only_on_cpu(scratch_memory_north, 1:nfreqj(j), 1, j), nlons[j], 1)
                         for j in 1:nlat_half]

    # BATCHED (K>1, 2-D) plans — a `Dict{Int, Vector{P2}}` keyed by K with CONCRETE plan type P2.
    # The Dict value type is fixed up-front from a throwaway 2-column plan (a plan's Julia type depends
    # only on eltype/ndims/region, not on the array size), so the Dict stays concretely typed even when
    # empty (nlayers==1). 
    NF_real = eltype(fake_grid_data.data)
    tmp_real = similar(fake_grid_data.data, NF_real, (nlons[1], 2))
    tmp_complex = similar(scratch_memory_north, Complex{NF_real}, (nfreqj(1), 2))
    RVec = Vector{typeof(AbstractFFTs.plan_rfft(tmp_real, 1))}
    BVec = Vector{typeof(AbstractFFTs.plan_brfft(tmp_complex, nlons[1], 1))}
    rfft_plans_batched = Dict{Int, RVec}()
    brfft_plans_batched = Dict{Int, BVec}()
    for K in planned_K
        K > 1 || continue                          # K=1 handled by the serial plans above
        rfft_plans_batched[K] = [AbstractFFTs.plan_rfft(view_only_on_cpu(fake_grid_data.data, rings[j], 1:K), 1)
                                 for j in 1:nlat_half]
        brfft_plans_batched[K] = [AbstractFFTs.plan_brfft(view_only_on_cpu(scratch_memory_north, 1:nfreqj(j), 1:K, j), nlons[j], 1)
                                  for j in 1:nlat_half]
    end

    return rfft_plan_serial, brfft_plan_serial, rfft_plans_batched, brfft_plans_batched
end
