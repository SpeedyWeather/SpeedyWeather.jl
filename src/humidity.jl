# These numbers are hard-coded in the original fortran function
# TODO(alistair): make this better - maybe move them into constants.jl?
const e0 = 6.108e-3  # TODO(alistair): Check the precise meaning of the fortran literal: 6.108e-3_p
const C₁ = 17.269
const C₂ = 21.875
const T₀ = 273.16
const T₁ = 35.86
const T₂ = 7.66

"""
Compute the saturation specific humidity for a single atmospheric level using the
    Clausius-Clapeyron relation.
"""
function get_saturation_specific_humidity(
    T::Array{NF,2},  # Absolute temperature for a single atmospheric level
    P::Array{NF,2},  # Surface pressure
    σ::NF,           # σ at the current atmospheric level
    M::Model,
) where {NF<:AbstractFloat}
    @unpack nlon, nlat = M.geometry

    Qsat = zeros(nlon, nlat)  # Saturation specific humidity

    @inbounds for j = 1:nlat, i = 1:nlon
        if T[i, j] > T₀
            Qsat[i, j] = e0 * exp(C₁ * (T[i, j] - T₀) / (T[i, j] - T₁))
        else
            Qsat[i, j] = e0 * exp(C₂ * (T[i, j] - T₀) / (T[i, j] - T₂))
        end
    end

    if σ <= 0
        @. Qsat = 622.0 * Qsat / (P[1, 1] - 0.378 * Qsat)
    else
        @. Qsat = 622.0 * Qsat / (σ * P - 0.378 * Qsat)
    end

    return Qsat
end

"""
Compute the saturation specific humidity using the Clausius-Clapeyron relation.
"""
function get_saturation_specific_humidity(
    T::Array{NF,3},  # Absolute temperature
    P::Array{NF,2},  # Surface pressure
    M::Model,
) where {NF<:AbstractFloat}
    @unpack nlon, nlat, nlev, σ_levels_full = M.geometry

    Qsat = zeros(nlon, nlat, nlev)  # Saturation specific humidity

    @inbounds for k = 1:nlev
        σ = σ_levels_full[k]
        Qsat[:, :, k] = get_saturation_specific_humidity(T[:, :, k], P, σ, M)
    end

    return Qsat
end

"""
Convert humidity into relative humidity.
"""
function get_relative_humidity(
    Q::Array{NF,3},  # Specific humidity
    T::Array{NF,3},  # Absolute temperature
    P::Array{NF,2},  # Surface pressure
    M::Model,
) where {NF<:AbstractFloat}
    Qsat = get_saturation_specific_humidity(T, P, M)
    return Q ./ Qsat
end
