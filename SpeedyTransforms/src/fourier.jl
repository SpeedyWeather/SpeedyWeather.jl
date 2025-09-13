# function barrier for batched or serial transforms as FFTW plans cannot be reused for fewer vertical layers
function _fourier!(f_north, f_south, field::AbstractField, S::SpectralTransform)
    _fourier! = if size(field, 2) == S.nlayers > 1
        _fourier_batched!
    else
        _fourier_serial!
    end
    return _fourier!(f_north, f_south, field, S)
end

# function barrier for batched or serial transforms as FFTW plans cannot be reused for fewer vertical layers
function _fourier!(field::AbstractField, f_north, f_south, S::SpectralTransform)
    _fourier! = if size(field, 2) == S.nlayers > 1
        _fourier_batched!
    else
        _fourier_serial!
    end
    return _fourier!(field, f_north, f_south, S)
end

"""$(TYPEDSIGNATURES)
Apply FFT plan to field. Due to different limitations of FFTW and CUDA FFT plans,
this function is dispatched on the array type and indexing.

With FFTW currently strided inputs are possible, there's low-level support for 
strided outputs, but it's not really implementend to be used at a high-level, so we 
use broadcasting instead of mul!
"""
function _apply_fft_plan!(f_out, nfreq, layers, j, plan, field::Field{NF, N, <:Array}, ilons, k_grid=Colon()) where {NF, N} 
    view(f_out, 1:nfreq, layers, j) .= plan * view(field.data, ilons, k_grid)
end
function _apply_fft_plan!(field_out, ilons, k, plan, g_in::Array, nfreq, layers, j) 
    view(field_out, ilons, k) .= plan * view(g_in, 1:nfreq, layers, j)
end

"""$(TYPEDSIGNATURES)
(Forward) FFT, applied in zonal direction of `field` provided. 
"""
function _apply_batched_fft!(
    f_out::AbstractArray{<:Complex, 3},
    field::AbstractField,
    S::SpectralTransform, 
    j::Int,
    nfreq::Int,
    ilons::UnitRange{Int};
    not_equator::Bool = true
)
    rfft_plan = S.rfft_plans[j]     # FFT planned wrt nlon on ring
    nlayers = size(field, 2)        # number of vertical layers

    if not_equator
        _apply_fft_plan!(f_out, nfreq, 1:nlayers, j, rfft_plan, field, ilons)
    else
        fill!(f_out[1:nfreq, 1:nlayers, j], 0)
    end
end

"""$(TYPEDSIGNATURES)
(Inverse) FFT, applied in zonal direction of `field` provided.
"""
function _apply_batched_fft!(
    field::AbstractField,
    g_in::AbstractArray{<:Complex, 3},
    S::SpectralTransform,
    j::Int,
    nlon::Int,
    ilons::UnitRange{Int};
    not_equator::Bool = true
)
    brfft_plan = S.brfft_plans[j]   # FFT planned wrt nlon on ring
    nlayers = size(field, 2)        # number of vertical layers
    nfreq = nlon÷2 + 1              # linear max Fourier frequency wrt to nlon

    if not_equator
        _apply_fft_plan!(field.data, ilons, Colon(), brfft_plan, g_in, nfreq, 1:nlayers, j)
    end
end

"""$(TYPEDSIGNATURES)
(Forward) FFT, applied in vertical direction of `field` provided. 
"""
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
    rfft_plan = S.rfft_plans_1D[j]     # FFT planned wrt nlon on ring
    k_grid = eachlayer(field)[k]        # vertical layer index

    if not_equator
        _apply_fft_plan!(f_out, nfreq, k, j, rfft_plan, field, ilons, k_grid)
    else
        fill!(f_out[1:nfreq, k, j], 0)
    end
end

"""$(TYPEDSIGNATURES)
(Inverse) FFT, applied in vertical direction of `field` provided.
"""
function _apply_serial_fft!(
    field::AbstractField,
    g_in::AbstractArray{<:Complex, 3},
    S::SpectralTransform,
    j::Int,
    k::Int,
    nfreq::Int,
    ilons::UnitRange{Int};
    not_equator::Bool = true
)
    brfft_plan = S.brfft_plans_1D[j]   # FFT planned wrt nlon on ring
    k_grid = eachlayer(field)[k]     # vertical layer index

    if not_equator
        _apply_fft_plan!(field.data, ilons, k_grid, brfft_plan, g_in, nfreq, k_grid, j)
    end
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
    S::SpectralTransform;                   # precomputed transform
)
    (; nlat, nlons) = S                     # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers == S.nlayers || throw(DimensionMismatch(S, field))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    rings = eachring(field)                 # precomputed ring indices
    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # corresponding southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring
        nfreq = nlon÷2 + 1                  # linear max Fourier frequency wrt to nlon
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        ilons = rings[j_north]              # in-ring indices northern ring
        # FOURIER TRANSFORM in zonal direction, northern latitude
        _apply_batched_fft!(f_north, field, S, j, nfreq, ilons)

        # and southern latitude if not on Equator
        ilons = rings[j_south]              # in-ring indices southern ring
        _apply_batched_fft!(f_south, field, S, j, nfreq, ilons; not_equator=not_equator)
       
    end
end

"""$(TYPEDSIGNATURES)
(Forward) Fast Fourier transform (grid to spectral) in zonal direction of `grids`,
stored in scratch memories `f_north`, `f_south` to be passed on to the Legendre transform.
Serial version that does not require the number of vertical layers to be the same as precomputed in `S`.
Not to be called directly, use `transform!` instead."""
function _fourier_serial!(                  # GRID TO SPECTRAL
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
    f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    field::AbstractField,                   # gridded input
    S::SpectralTransform;                   # precomputed transform
)
    (; nlat, nlons) = S                     # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    rings = eachring(field)                     # precomputed ring indices
    @inbounds for (k, k_grid) in zip(1:nlayers, eachlayer(field))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j = j_north                         # symmetric index / ring-away from pole index
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j]                     # number of longitudes on this ring (north or south)
            nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            ilons = rings[j_north]              # in-ring indices northern ring
            
            # Apply FFT in the northern latitudes
            _apply_serial_fft!(f_north, field, S, j, k, nfreq, ilons)
            
            # southern latitude, don't call redundant 2nd fft if ring is on equator 
            ilons = rings[j_south]                      # in-ring indices southern ring
            _apply_serial_fft!(f_south, field, S, j, k, nfreq, ilons; not_equator=not_equator)
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
    S::SpectralTransform;                   # precomputed transform
)
    (; nlat, nlons) = S                     # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers == S.nlayers || throw(DimensionMismatch(S, field))     # otherwise FFTW complains
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    rings = eachring(field)                 # precomputed ring indices
    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring (north or south)
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        ilons = rings[j_north]              # in-ring indices northern ring

        # northern latitude
        _apply_batched_fft!(field, g_north, S, j, nlon, ilons)

        # southern latitude, don't call redundant 2nd FFT if ring is on equator 
        ilons = rings[j_south]              # in-ring indices southern ring
        _apply_batched_fft!(field, g_south, S, j, nlon, ilons; not_equator=not_equator)

    end
end

"""$(TYPEDSIGNATURES)
(Inverse) Fast Fourier transform (spectral to grid) of Legendre-transformed inputs `g_north` and `g_south`
to be stored in `grids`. Serial version that does not require the number of vertical layers to be the same
as precomputed in `S`. Not to be called directly, use `transform!` instead."""
function _fourier_serial!(                  # SPECTRAL TO GRID
    field::AbstractField,                   # gridded output
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
    g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    S::SpectralTransform;                   # precomputed transform
)
    (; nlat, nlons) = S                     # dimensions   
    (; nlat_half) = S.grid       
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))     # otherwise FFTW complains
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    rings = eachring(field)                     # precomputed ring indices
    @inbounds for (k, k_grid) in zip(1:nlayers, eachlayer(field))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j = j_north                         # symmetric index / ring-away from pole index
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j]                     # number of longitudes on this ring (north or south)
            nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            # Apply FFT in the northern latitudes
            ilons = rings[j_north]              # in-ring indices northern ring
            _apply_serial_fft!(field, g_north, S, j, k, nfreq, ilons)

            # southern latitude, don't call redundant 2nd fft if ring is on equator
            ilons = rings[j_south]              # in-ring indices southern ring
            _apply_serial_fft!(field, g_south, S, j, k, nfreq, ilons; not_equator=not_equator)
        end
    end
end

which_FFT_package(::Type{<:AbstractArray{NF}}) where NF<:AbstractFloat = NF <: Union{Float32, Float64} ? FFTW : GenericFFT

"""$(TYPEDSIGNATURES)
Util function to generate FFT plans based on the array type of the fake Grid 
data provided. Uses views, which is less allocate-y than indexing but breaks 
when using CuArrays (see CUDA extension for alternative implementation for 
CuArrays)."""
function plan_FFTs!(
    rfft_plans::Vector{AbstractFFTs.Plan},
    brfft_plans::Vector{AbstractFFTs.Plan},
    rfft_plans_1D::Vector{AbstractFFTs.Plan},
    brfft_plans_1D::Vector{AbstractFFTs.Plan},
    fake_grid_data::AbstractField{NF, N, <:AbstractArray{NF}},
    scratch_memory_north::AbstractArray{Complex{NF}},
    rings,
    nlons::Vector{<:Int}
) where {NF<:AbstractFloat, N}
    # Determine which FFT package to use (currently either FFTW or GenericFFT)
    FFT_package = which_FFT_package(Array{NF})

    # For each ring generate an FFT plan (for all layers and for a single layer)
    for (j, nlon) in enumerate(nlons)
        real_matrix_input = view(fake_grid_data.data, rings[j], :)
        complex_matrix_input = view(scratch_memory_north, 1:nlon÷2 + 1, :, j)
        real_vector_input = view(fake_grid_data.data, rings[j], 1)
        complex_vector_input = view(scratch_memory_north, 1:nlon÷2 + 1, 1, j)

        rfft_plans[j] = FFT_package.plan_rfft(real_matrix_input, 1)
        brfft_plans[j] = FFT_package.plan_brfft(complex_matrix_input, nlon, 1)
        rfft_plans_1D[j] = FFT_package.plan_rfft(real_vector_input, 1)
        brfft_plans_1D[j] = FFT_package.plan_brfft(complex_vector_input, nlon, 1)
    end

    return rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D
end