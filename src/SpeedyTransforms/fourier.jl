"""$(TYPEDSIGNATURES)
Inverse fast Fourier transform of Legendre-transformed inputs `g_north` and `g_south` to be stored in `grids`.
Not to be called directly, use `transform!` instead."""
function _fourier!(
    grids::AbstractGridArray,               # gridded output
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
    g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat, nlons, nlat_half, nlayers) = S # dimensions
    (; brfft_plans) = S                     # pre-planned transforms

    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, nlayers, nlat_half) || throw(DimensionMismatch(S, grids))

    rings = eachring(grids)                 # precomputed ring indices

    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index
        j_south = nlat - j_north + 1        # southern latitude index
        nlon = nlons[j]                     # number of longitudes on this ring (north or south)
        nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        brfft_plan = brfft_plans[j]         # FFT planned wrt nlon on ring
        ilons = rings[j_north]              # in-ring indices northern ring

        # PERFORM FFT, inverse complex to real, hence brfft
        # FFTW is in-place writing into `out` via `mul`
        # FFTW requires stride-1 output, hence view on scratch memory
        out = reshape(view(S.scratch_memory_grid, 1:nlon*nlayers), (nlon, nlayers))
        LinearAlgebra.mul!(out, brfft_plan, view(g_north, 1:nfreq, :, j))
        grids[ilons, :] = out

        # southern latitude, don't call redundant 2nd FFT if ring is on equator 
        ilons = rings[j_south]              # in-ring indices southern ring
        if not_equator
            LinearAlgebra.mul!(out, brfft_plan, view(g_south, 1:nfreq, :, j))
            grids[ilons, :] = out
        end
    end
end