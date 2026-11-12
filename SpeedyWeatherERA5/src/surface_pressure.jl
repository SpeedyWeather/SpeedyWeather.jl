# Surface pressure for the model. ERA5's ps belongs to ERA5's orography; the model has its own (smoothed) orography.
#   raw_pressure                      ln ps truncated to the model resolution, no adjustment (the simple option)
#   orography_consistent_pressure     hydrostatic shift ln ps_model = ln ps_ERA5 + g (h_ERA5 - h_model)/(R Ts), done as
#                                     X = ln ps + g h_ERA5/(R Ts) at native resolution, truncate X, then ln ps = X - g h_model/(R Ts)
#                                     (truncating the smooth X rather than ln ps avoids Gibbs ringing at orography)

"""Temperature at the surface pressure `ps`: lowest above-ground ERA5 level extrapolated with the standard lapse rate."""
function surface_temperature(T, lev, ps; R = 287.0f0, g = 9.81f0)
    Ts = similar(ps)
    κΓ = R * Γ_STD / g
    Threads.@threads for j in axes(ps, 2)
        for i in axes(ps, 1)
            nv = max(searchsortedlast(lev, ps[i, j]), 1)
            Ts[i, j] = T[i, j, nv] * (ps[i, j] / lev[nv])^κΓ
        end
    end
    return Ts
end

"""Naive option: ERA5's `ln ps` truncated to the model resolution and put on the model orography unchanged."""
raw_pressure(on_model, ps) = on_model(log.(ps))

"""
`ln ps` [Pa] on the model grid consistent with the MODEL's orography: reduce ERA5's `ps` with ERA5's orography `zs` (geopotential),
truncate to the model resolution (`on_model`), and rebuild with `model.orography`. `Ts` is the surface temperature, the lowest
above-ground ERA5 level extrapolated with the standard lapse rate (floored at 200 K after truncation)."""
function orography_consistent_pressure(model, on_model, T, lev, ps, zs; R, g)
    Ts = surface_temperature(T, lev, ps; R, g)
    X = log.(ps) .+ g .* (zs ./ g) ./ (R .* Ts)                 # smooth: orography signal removed before truncating
    X63, Ts63 = on_model(X), max.(on_model(Ts), 200.0f0)
    return X63 .- g .* model.orography.orography ./ (R .* Ts63)
end
