"""$(TYPEDSIGNATURES)
Fast Fourier transform in zonal direction of `grids`.""" 
function _fourier!(                         # GRID TO SPECTRAL
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
    f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    grids::AbstractGridArray,               # gridded output
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat, nlons, nlat_half) = S          # dimensions
    (; rfft_plans) = S                      # pre-planned transforms
    nlayers = size(grids, 2)                # number of vertical layers

    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
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
        grid_in = reshape(view(S.scratch_memory_grid, 1:nlon*nlayers), (nlon, nlayers))
        grid_in[:,:] = view(grids.data, ilons, :)
        out = view(f_north, 1:nfreq, 1:nlayers, j)
        LinearAlgebra.mul!(out, rfft_plan, grid_in)   # Northern latitude

        # and southern latitude if not on Equator (
        ilons = rings[j_south]              # in-ring indices southern ring
        grid_in[:,:] = view(grids.data, ilons, :)
        out = view(f_south, 1:nfreq, 1:nlayers, j)
        if not_equator                      # skip FFT, redundant because north did that latitude already
            LinearAlgebra.mul!(out, rfft_plan, grid_in)
        else
            out .= 0    # TODO necessary?
        end
    end
end

"""$(TYPEDSIGNATURES)
Inverse fast Fourier transform of Legendre-transformed inputs `g_north` and `g_south` to be stored in `grids`.
Not to be called directly, use `transform!` instead."""
function _fourier!(                         # SPECTRAL TO GRID
    grids::AbstractGridArray,               # gridded output
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
    g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat, nlons, nlat_half) = S          # dimensions
    (; brfft_plans) = S                     # pre-planned transforms
    nlayers = size(grids, 2)                # number of vertical layers

    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, grids))

    rings = eachring(grids)                 # precomputed ring indices

    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring (north or south)
        nfreq = nlon÷2 + 1                  # linear max Fourier frequency wrt to nlon
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        brfft_plan = brfft_plans[j]         # FFT planned wrt nlon on ring
        ilons = rings[j_north]              # in-ring indices northern ring

        # PERFORM FFT, inverse complex to real, hence brfft
        # FFTW is in-place writing into `out` via `mul`
        # FFTW requires stride-1 output, hence view on scratch memory
        out = reshape(view(S.scratch_memory_grid, 1:nlon*nlayers), (nlon, nlayers))
        LinearAlgebra.mul!(out, brfft_plan, view(g_north, 1:nfreq, 1:nlayers, j))
        grids[ilons, :] = out

        # southern latitude, don't call redundant 2nd FFT if ring is on equator 
        ilons = rings[j_south]              # in-ring indices southern ring
        if not_equator
            LinearAlgebra.mul!(out, brfft_plan, view(g_south, 1:nfreq, 1:nlayers, j))
            grids[ilons, :] = out
        end
    end
end