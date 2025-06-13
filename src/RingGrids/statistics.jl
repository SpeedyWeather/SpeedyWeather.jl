"""$(TYPEDSIGNATURES)
Zonal mean of `grid`, i.e. along its latitude rings."""
function zonal_mean(field::AbstractField)
    # extend to any non-horizontal dimensions of grid (e.g. vertical or time)
    ks = size(field)[2:end]

    # determine type T after division with integer (happening in mean)
    T = Base.promote_op(/, eltype(field), Int64)
    m = zeros(T, get_nlat(field), ks...)

    rings = eachring(field.grid)
    for k in eachlayer(field)
        for (j, ring) in enumerate(rings)
            m[j, k] = mean(field[ring, k])
        end
    end

    return m
end