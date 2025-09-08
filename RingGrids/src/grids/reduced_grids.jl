# only defined for full grids, return empty arrays here
get_lon(::Type{<:AbstractReducedGrid}, nlat_half::Integer) = Float64[]
get_lond(::Type{<:AbstractReducedGrid}, nlat_half::Integer) = Float64[]

# all reduced grids have their maximum number of longitude points around the equator, i.e. j = nlat_half
get_nlon_max(Grid::Type{<:AbstractReducedGrid}, nlat_half::Integer) = get_nlon_per_ring(Grid, nlat_half, nlat_half)

"""$(TYPEDSIGNATURES)
Return vectors of `londs`, `latds`, of longitudes (degrees, 0-360˚E)
and latitudes (degrees, -90˚ to 90˚N) for reduced grids for all
grid points in order of the running index `ij`."""
function get_londlatds(Grid::Type{<:AbstractReducedGrid}, nlat_half::Integer)
    
    latd = get_latd(Grid, nlat_half)
    nlat = get_nlat(Grid, nlat_half)
    
    npoints = get_npoints(Grid, nlat_half)
    londs = zeros(npoints)      # preallocate arrays
    latds = zeros(npoints)

    ij = 1                      # running index
    for j in 1:nlat             # populate arrays ring by ring
        lond = get_lond_per_ring(Grid, nlat_half, j)
        nlon = length(lond)

        latds[ij:ij+nlon-1] .= latd[j]
        londs[ij:ij+nlon-1] .= lond

        ij += nlon
    end

    return londs, latds
end