"""$(TYPEDSIGNATURES)
Zonal mean of `grid`, i.e. along its latitude rings."""
function zonal_mean(grid::AbstractGridArray)
    # extend to any non-horizontal dimensions of grid (e.g. vertical or time)
    ks = size(grid)[2:end]
    m = zeros(eltype(grid), RingGrids.get_nlat(grid), ks...)

    rings = RingGrids.eachring(grid)
    for k in eachgrid(grid)
        for (j, ring) in enumerate(rings)
            m[j, k] = mean(grid[ring, k])
        end
    end

    return m
end