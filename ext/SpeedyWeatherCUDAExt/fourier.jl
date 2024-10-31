
# Override FFT package deciding function
SpeedyTransforms.which_FFT_package(::Type{<:CuArray{<:AbstractFloat}}) = CUFFT 

"""$(TYPEDSIGNATURES)
Util function to generate FFT plans based on the array type of the fake Grid 
data provided. Uses indexing as we seemingly can't use views for FFT planning 
with CUFFT."""
function SpeedyTransforms.plan_FFTs!(
    rfft_plans::Vector{AbstractFFTs.Plan},
    brfft_plans::Vector{AbstractFFTs.Plan},
    rfft_plans_1D::Vector{AbstractFFTs.Plan},
    brfft_plans_1D::Vector{AbstractFFTs.Plan},
    fake_grid_data::AbstractGridArray{NF, N, <:CuArray{NF}},
    scratch_memory_north::CuArray{Complex{NF}},
    rings::AbstractArray,
    nlons::Vector{<:Int}
) where {NF<:AbstractFloat, N}
    # Determine which FFT package to use (currently either FFTW or GenericFFT)
    FFT_package = SpeedyTransforms.which_FFT_package(CuArray{NF})

    # For each ring generate an FFT plan (for all layers and for a single layer)
    for (j, nlon) in enumerate(nlons)
        real_matrix_input = fake_grid_data.data[rings[j], :]
        complex_matrix_input = scratch_memory_north[1:nlon÷2 + 1, :, j]
        real_vector_input = fake_grid_data.data[rings[j], 1]
        complex_vector_input = scratch_memory_north[1:nlon÷2 + 1, 1, j]

        rfft_plans[j] = FFT_package.plan_rfft(real_matrix_input, 1)
        brfft_plans[j] = FFT_package.plan_brfft(complex_matrix_input, nlon, 1)
        rfft_plans_1D[j] = FFT_package.plan_rfft(real_vector_input, 1)
        brfft_plans_1D[j] = FFT_package.plan_brfft(complex_vector_input, nlon, 1)
    end

    return rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D
end


function SpeedyTransforms._fourier_batched!(                     # GRID TO SPECTRAL
    f_north::CuArray{<:Complex, 3},             # Fourier-transformed output
    f_south::CuArray{<:Complex, 3},             # and for southern latitudes
    grids::AbstractGridArray{NF, N, <:CuArray}, # gridded input
    S::SpectralTransform,                       # precomputed transform
) where {NF<:AbstractFloat, N}
    (; nlat, nlons, nlat_half) = S          # dimensions
    (; rfft_plans) = S                      # pre-planned transforms
    nlayers = size(grids, 2)

    @boundscheck SpeedyTransforms.ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck nlayers == S.nlayers || throw(DimensionMismatch(S, grids))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, grids))

    rings = eachring(grids)                 # precomputed ring indices
    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # corresponding southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring
        nfreq = nlon÷2 + 1                  # linear max Fourier frequency wrt to nlon
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        rfft_plan = rfft_plans[j]           # FFT planned wrt nlon on ring
        ilons = rings[j_north]              # in-ring indices northern ring
        
        # FOURIER TRANSFORM in zonal direction, northern latitude
        view(f_north, 1:nfreq, 1:nlayers, j) .= rfft_plan * grids.data[ilons, :]

        # This is faster but doesn't seem to work for greater than 1D, seems to 
        # hang.
        # view(f_north, 1:nfreq, 1:nlayers, j) .= rfft_plan * view(grids.data, ilons, :)

        # and southern latitude if not on Equator
        ilons = rings[j_south]              # in-ring indices southern ring
        if not_equator                      # skip FFT, redundant because north did that latitude already
            view(f_south, 1:nfreq, 1:nlayers, j) .= rfft_plan * grids.data[ilons, :]
        else
            fill!(f_south[1:nfreq, 1:nlayers, j], 0)
        end
    end

end

function SpeedyTransforms._fourier_serial!(                      # GRID TO SPECTRAL
    f_north::CuArray{<:Complex, 3},             # Fourier-transformed output
    f_south::CuArray{<:Complex, 3},             # and for southern latitudes
    grids::AbstractGridArray{NF, N, <:CuArray}, # gridded input
    S::SpectralTransform,                       # precomputed transform
) where {NF<:AbstractFloat, N}
    (; nlat, nlons, nlat_half) = S          # dimensions
    rfft_plans = S.rfft_plans_1D            # pre-planned transforms
    nlayers = size(grids, 2)                # number of vertical layers

    @boundscheck SpeedyTransforms.ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, grids))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, grids))

    rings = eachring(grids)                     # precomputed ring indices
    @inbounds for (k, k_grid) in zip(1:nlayers, eachgrid(grids))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j = j_north                         # symmetric index / ring-away from pole index
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j]                     # number of longitudes on this ring (north or south)
            nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            # FOURIER TRANSFORM in zonal direction
            rfft_plan = rfft_plans[j]           # FFT planned wrt nlon on ring
            ilons = rings[j_north]              # in-ring indices northern ring
            view(f_north, 1:nfreq, k, j) .= rfft_plan * view(grids.data, ilons, k_grid)
            
            # southern latitude, don't call redundant 2nd fft if ring is on equator 
            ilons = rings[j_south]                      # in-ring indices southern ring
            if not_equator 
                view(f_south, 1:nfreq, k, j) .= rfft_plan * view(grids.data, ilons, k_grid) # perform FFT
            else
                fill!(view(f_south, 1:nfreq, k, j), 0)
            end
        end
    end
end
