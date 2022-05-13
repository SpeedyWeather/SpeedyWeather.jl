# These numbers are hard-coded in the original fortran function
# TODO(alistair): make this better - maybe move them into constants.jl?
const e₀ = 6.108e-3  # TODO(alistair): Check the precise meaning of the fortran literal: 6.108e-3_p
const C₁ = 17.269
const C₂ = 21.875
const T₀ = 273.16
const T₁ = 35.86
const T₂ = 7.66

"""
Compute the saturation specific humidity for a single atmospheric level.
"""
function get_saturation_specific_humidity(
    T::Array{NF,2},  # Absolute temperature for a single atmospheric level
    p::Array{NF,2},  # Normalised surface pressure
    σ::NF,           # σ of the current atmospheric level
    M::Model,
) where {NF<:AbstractFloat}
    @unpack nlon, nlat = M.geometry

    Qsat = zeros(nlon, nlat)  # Saturation specific humidity

    @inbounds for j = 1:nlat, i = 1:nlon
        if T[i, j] > T₀
            Qsat[i, j] = e₀ * exp(C₁ * (T[i, j] - T₀) / (T[i, j] - T₁))
        else
            Qsat[i, j] = e₀ * exp(C₂ * (T[i, j] - T₀) / (T[i, j] - T₂))
        end
    end

    if σ <= 0
        @. Qsat = 622.0 * Qsat / (p[1, 1] - 0.378 * Qsat)
    else
        @. Qsat = 622.0 * Qsat / (σ * p - 0.378 * Qsat)
    end

    return Qsat
end

"""
Compute the saturation specific humidity.
"""
function get_saturation_specific_humidity(
    T::Array{NF,3},  # Absolute temperature
    p::Array{NF,2},  # Normalised surface pressure
    M::Model,
) where {NF<:AbstractFloat}
    @unpack nlon, nlat, nlev, σ_levels_full = M.geometry

    Qsat = zeros(nlon, nlat, nlev)  # Saturation specific humidity

    @inbounds for k = 1:nlev
        σ = σ_levels_full[k]
        Qsat[:, :, k] = get_saturation_specific_humidity(T[:, :, k], p, σ, M)
    end

    return Qsat
end

"""
Convert specific humidity into relative humidity.
"""
function get_relative_humidity(
    Q::Array{NF,3},  # Specific humidity
    T::Array{NF,3},  # Absolute temperature
    p::Array{NF,2},  # Normalised surface pressure
    M::Model,
) where {NF<:AbstractFloat}
    Qsat = get_saturation_specific_humidity(T, p, M)
    return Q ./ Qsat
end
