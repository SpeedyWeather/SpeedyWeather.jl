"""
Compute spectral geopotential from spectral temperature `Tabs`
and spectral topography `ϕ0`.
"""
function geopotential!( ϕ::Array{T,3},      # geopotential
                        ϕ0::Array{T,2},     # geop
                        Tabs::Array{T,3},   # absolute Temperature
                        geometry::Geometry{T}) where {T<:AbstractFloat}

    mx,nx,nlev = size(ϕ)

    @boundscheck size(ϕ) == size(Tabs) || throw(BoundsError())
    @boundscheck (mx,nx) == size(ϕ0)   || throw(BoundsError())

    @unpack xgeop1, xgeop2, lapserate_correction = geometry

    # Bottom layer (integration over half a layer) is last index
    ϕ[:,:,end] = ϕ0 + xgeop1[end]*Tabs[:,:,end]

    # Other layers (integrate two half-layers from bottom to top)
    for k in nlev-1:-1:1
        ϕ[:,:,k] = ϕ[:,:,k+1] + xgeop2[k+1]*Tabs[:,:,k+1] + xgeop1[k]*Tabs[:,:,k]
    end

    # Lapse-rate correction in the free troposphere
    # TODO only for spectral coefficients 1,: ?
    for k in 2:nlev-1
        ϕ[1,:,k] = ϕ[1,:,k] +
                    lapserate_correction[k-1]*(Tabs[1,:,k+1] - Tabs[1,:,k-1])
    end
end
