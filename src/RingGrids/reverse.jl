Base._reverse!(field::AbstractField, sym::Symbol) = Base._reverse!(field, Val(sym))

function Base._reverse!(field::AbstractField, ::Val{:lat})

    rings = eachring(field)
    nrings = length(rings)
    rings_north = view(rings, 1:nrings ÷ 2)

    @inbounds for k in eachfield(field)
        for (j_north, ring) in enumerate(rings_north)
            j_south = nrings - j_north + 1
            for (i, ij_north) in enumerate(ring)
                ij_south = rings[j_south][i]

                # read out northern and southern values
                north = field[ij_north, k]
                south = field[ij_south, k]

                # and swap them
                field[ij_north, k] = south
                field[ij_south, k] = north
            end
        end
    end

    return field
end

function Base._reverse!(field::AbstractField, ::Val{:lon})
    for k in eachfield(field)
        for ring in eachring(field)
            field_ring = view(field, ring, k)
            reverse!(field_ring)
        end
    end
    return field
end

# also allow long names longitude, latitude
Base._reverse!(field::AbstractField, ::Val{:latitude}) =  Base._reverse!(field, Val{:lat}())
Base._reverse!(field::AbstractField, ::Val{:longitude}) = Base._reverse!(field, Val{:lon}())