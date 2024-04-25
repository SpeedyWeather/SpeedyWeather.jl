# functions to scale the latitude of any grid, in-place
scale_coslat!(  grid::AbstractGridArray) = _scale_coslat!(grid, power=1)
scale_coslat²!( grid::AbstractGridArray) = _scale_coslat!(grid, power=2)
scale_coslat⁻¹!(grid::AbstractGridArray) = _scale_coslat!(grid, power=-1)
scale_coslat⁻²!(grid::AbstractGridArray) = _scale_coslat!(grid, power=-2)

# and via a deepcopy
scale_coslat(  grid::AbstractGridArray) = scale_coslat!(  deepcopy(grid))
scale_coslat²( grid::AbstractGridArray) = scale_coslat²!( deepcopy(grid))
scale_coslat⁻¹(grid::AbstractGridArray) = scale_coslat⁻¹!(deepcopy(grid))
scale_coslat⁻²(grid::AbstractGridArray) = scale_coslat⁻²!(deepcopy(grid))

function _scale_coslat!(grid::Grid; power=1) where {Grid<:AbstractGridArray}
    lat = get_lat(Grid, grid.nlat_half)
    T = eltype(grid)
    coslat = @. convert(T, cos(lat)^power)
    return _scale_lat!(grid, coslat)
end

"""
$(TYPEDSIGNATURES)
Generic latitude scaling applied to `A` in-place with latitude-like vector `v`."""
function _scale_lat!(grid::AbstractGridArray{T}, v::AbstractVector) where T
    @boundscheck get_nlat(grid) == length(v) || throw(BoundsError)
    
    @inbounds for k in eachgrid(grid)
        for (j, ring) in enumerate(eachring(grid))
            vj = convert(T, v[j])
            for ij in ring
                grid[ij, k] *= vj
            end
        end
    end

    return grid
end