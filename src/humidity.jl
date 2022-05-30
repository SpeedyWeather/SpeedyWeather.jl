"""
Compute the saturation vapour pressure.
"""
function get_saturation_vapour_pressure!(
    sat_vap_pressure ::Array{NF,3},
    temp_grid        ::Array{NF,3},
    model            ::ModelSetup,
    ) where {NF<:AbstractFloat}
    @unpack nlon, nlat, nlev = model.geospectral.geometry
    @unpack e₀, C₁, C₂, T₀, T₁, T₂ = model.parameters.humidity_coefs

    for k = 1:nlev, j = 1:nlat, i = 1:nlon
        if temp_grid[i, j, k] > T₀
            # Saturation vapour pressure over water
            sat_vap_pressure[i, j, k] = e₀ * exp(C₁ * (temp_grid[i, j, k] - T₀) / (temp_grid[i, j, k] - T₁))
        else
            # Saturation vapour pressure over ice
            sat_vap_pressure[i, j, k] = e₀ * exp(C₂ * (temp_grid[i, j, k] - T₀) / (temp_grid[i, j, k] - T₂))
        end
    end

    return nothing
end

"""
Compute the saturation specific humidity.
"""
function get_saturation_specific_humidity!(
    sat_spec_humidity ::Array{NF,3},
    sat_vap_pressure  ::Array{NF,3},
    temp_grid         ::Array{NF,3},
    pres              ::Array{NF,2},
    model             ::ModelSetup,
) where {NF<:AbstractFloat}
    @unpack nlon, nlat, nlev, σ_levels_full = model.geospectral.geometry

    get_saturation_vapour_pressure!(sat_vap_pressure, temp_grid, model)

    for k = 1:nlev, j = 1:nlat, i = 1:nlon
        sat_spec_humidity[i, j, k] = 0.622 * sat_vap_pressure[i, j, k] / (pres[i, j] * σ_levels_full[k] - 0.378 * sat_vap_pressure[i, j, k])
    end

    return nothing
end
