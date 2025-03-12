Base._reverse!(grid::AbstractGridArray, sym::Symbol) = Base._reverse!(grid, Val(sym))
Base._reverse!(grid::AbstractGridArray, ::Val{Colon()}) = reverse!(grid)

function Base._reverse!(grid::AbstractGridArray, ::Val{:lat})

    rings = eachring(grid)
    nrings = length(rings)
    rings_north = view(rings, 1:nrings ÷ 2)

    @inbounds for k in eachgrid(grid)
        for (j_north, ring) in enumerate(rings_north)
            j_south = nrings - j_north + 1
            for (i, ij_north) in enumerate(ring)
                ij_south = rings[j_south][i]

                # read out northern and southern values
                north = grid[ij_north, k]
                south = grid[ij_south, k]

                # and swap them
                grid[ij_north, k] = south
                grid[ij_south, k] = north
            end
        end
    end

    return grid
end

function Base._reverse!(grid::AbstractGridArray, ::Val{:lon})
    for k in eachgrid(grid)
        for ring in eachring(grid)
            grid_ring = view(grid, ring, k)
            reverse!(grid_ring)
        end
    end
    return grid
end

# also allow long names longitude, latitude
Base._reverse!(grid::AbstractGridArray, ::Val{:latitude}) = Base._reverse!(grid, Val{:lat}())
Base._reverse!(grid::AbstractGridArray, ::Val{:longitude}) = Base._reverse!(grid, Val{:lon}())