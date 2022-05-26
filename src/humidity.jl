"""
Compute the saturation specific humidity.
"""
function get_saturation_specific_humidity!(
    T::Array{NF,3},  # Absolute temperature
    p::Array{NF,2},  # Normalised surface pressure
    M::ModelSetup{NF},
) where {NF<:AbstractFloat}
    @unpack nlon, nlat, nlev, σ_levels_full = M.geospectral.geometry
    @unpack humid_saturation = Diag.parametrization_variables
    @unpack e₀, C₁, C₂, T₀, T₁, T₂ = M.parameters.humidity_coefs

    for k = 1:nlev, j = 1:nlat, i = 1:nlon
        if T[i, j, k] > T₀
            Qsat[i, j, k] = e₀ * exp(C₁ * (T[i, j, k] - T₀) / (T[i, j, k] - T₁))
        else
            Qsat[i, j, k] = e₀ * exp(C₂ * (T[i, j, k] - T₀) / (T[i, j, k] - T₂))
        end
    end

    @. Qsat = 622.0 * Qsat / (σ_levels_full * p - 0.378 * Qsat)  # What does this do?

    return nothing
end

"""
Convert specific humidity into relative humidity.
"""
function get_relative_humidity(
    Q::Array{NF,3},  # Specific humidity
    T::Array{NF,3},  # Absolute temperature
    p::Array{NF,2},  # Normalised surface pressure
    M::ModelSetup{NF},
) where {NF<:AbstractFloat}
    Qsat = get_saturation_specific_humidity(T, p, M)
    return Q ./ Qsat
end
