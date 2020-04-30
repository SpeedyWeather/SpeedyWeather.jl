struct Boundaries{T<:AbstractFloat}
    ϕ0::Array{T,2}              # Surface geopotential (i.e. orography)
    ϕ0trunc::Array{T,2}         # Spectrally truncated surface geopotential
    land_sea_mask::Array{T,2}   # land-sea mask
    albedo::Array{T,2}          # Annual mean surface albedo
end

function Boundaries(T, geometry::Geometry, spectral_trans::SpectralTrans, g)
    # Read surface geopotential (i.e. orography)
    ϕ₀ = g*load_boundary_file("surface.nc", "orog")

    # Also store spectrally truncated surface geopotential
    ϕ₀ₛ = spectral_truncation(geometry, spectral_trans, ϕ₀)

    # Read land-sea mask
    land_sea_mask = load_boundary_file("surface.nc", "lsm")

    # Annual-mean surface albedo
    albedo = load_boundary_file("surface.nc", "alb")

    Boundaries(convert(Matrix{T}, ϕ₀), convert(Matrix{T}, ϕ₀ₛ), convert(Matrix{T}, land_sea_mask),
               convert(Matrix{T}, albedo))
end

# Compute a spectrally-filtered grid-point field.
function spectral_truncation(geometry::Geometry, spectral_trans::SpectralTrans, input::AbstractMatrix)
    input_sp = grid_to_spec(geometry, spectral_trans, input)

    for n in 1:geometry.nx
        for m in 1:geometry.mx
            N = m + n - 2
            if N > geometry.trunc
                input_sp[m,n] = Complex{typeof(input[1,1])}(0.0)
            end
        end
    end

    spec_to_grid(geometry, spectral_trans, input_sp)
end
