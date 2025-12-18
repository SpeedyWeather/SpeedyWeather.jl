module SpeedyTransformsCUDAExt

using CUDA
import CUDA: CUDA, CUFFT, CuArray, @captured
import AbstractFFTs
import LinearAlgebra
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays
import SpeedyTransforms.Architectures: GPU
import SpeedyTransforms: ismatching, _apply_batched_fft!, _apply_serial_fft!

"""$(TYPEDSIGNATURES)
(Forward) Fast Fourier transform (grid to spectral) in zonal direction of `field`,
stored in scratch memories `f_north`, `f_south` to be passed on to the Legendre transform.
Batched version that requires the number of vertical layers to be the same as precomputed in `S`.
Not to be called directly, use `transform!` instead."""
function SpeedyTransforms._fourier_batched!(                 # GRID TO SPECTRAL
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
    f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    field::AbstractField,                   # gridded input
    S::SpectralTransform{NF, <:GPU};                   # precomputed transform
) where NF
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers == S.nlayers || throw(DimensionMismatch(S, field))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    @inbounds @captured for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
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
(Forward) Fast Fourier transform (grid to spectral) in zonal direction of `field`,
stored in scratch memories `f_north`, `f_south` to be passed on to the Legendre transform.
Serial version that does not require the number of vertical layers to be the same as precomputed in `S`.
Not to be called directly, use `transform!` instead."""
function SpeedyTransforms._fourier_serial!(                  # GRID TO SPECTRAL
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed output
    f_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    field::AbstractField,                   # gridded input
    S::SpectralTransform{NF, <:GPU};                   # precomputed transform
) where NF
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need too match."
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    @inbounds @captured for (k, k_grid) in zip(1:nlayers, eachlayer(field))
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
function SpeedyTransforms._fourier_batched!(                 # SPECTRAL TO GRID
    field::AbstractField,                   # gridded output
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
    g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    S::SpectralTransform{NF, <:GPU};                   # precomputed transform
) where NF
    (; nlat, nlons, rings) = S              # dimensions
    (; nlat_half) = S.grid
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers == S.nlayers || throw(DimensionMismatch(S, field))     # otherwise FFTW complains
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    @inbounds @captured for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
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
to be stored in `field`. Serial version that does not require the number of vertical layers to be the same
as precomputed in `S`. Not to be called directly, use `transform!` instead."""
function SpeedyTransforms._fourier_serial!(                  # SPECTRAL TO GRID
    field::AbstractField,                   # gridded output
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed input
    g_south::AbstractArray{<:Complex, 3},   # and for southern latitudes
    S::SpectralTransform{NF, <:GPU};                   # precomputed transform
) where NF
    (; nlat, nlons, rings) = S              # dimensions   
    (; nlat_half) = S.grid       
    nlayers = size(field, 2)                # number of vertical layers

    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, field))     # otherwise FFTW complains
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, field))

    @inbounds @captured for (k, k_grid) in zip(1:nlayers, eachlayer(field))
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

end 
