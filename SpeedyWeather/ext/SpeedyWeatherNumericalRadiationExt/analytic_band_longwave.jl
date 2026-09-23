# NumericalRadiation's analytic-band longwave scheme (Williams 2026) as a SpeedyWeather
# longwave component: Radiation(spectral_grid; longwave = AnalyticBandLongwave(spectral_grid)).
# The scheme type is NumericalRadiation's own; this file adds the methods SpeedyWeather calls.

"""$(TYPEDSIGNATURES)
An `AnalyticBandLongwave` in the number format of the spectral grid; `kwargs` are the scheme's
parameters. CO₂ follows the model's `greenhouse_gases.co2` when present, else
$(DEFAULT_CO₂) ppm."""
AnalyticBandLongwave(SG::SpeedyWeather.SpectralGrid; kwargs...) = AnalyticBandLongwave{SG.NF}(; kwargs...)

# the standard longwave diagnostics, as declared for any SpeedyWeather.AbstractLongwave
SpeedyWeather.variables(::AnalyticBandLongwave{NF}, ::SpeedyWeather.AbstractModel) where NF =
    SpeedyWeather.variables(SpeedyWeather.UniformCooling{NF}())

SpeedyWeather.initialize!(::AnalyticBandLongwave, ::SpeedyWeather.PrimitiveEquation) = nothing

# Every constant comes from the SpeedyWeather model, so the radiation runs
# with the host's values; SpeedyWeather stores molar masses in g mol⁻¹.
@inline function speedy_physical_constants(model)
    NF = typeof(model.planet.gravity)
    (; planet, atmosphere) = model
    return PhysicalConstants{NF}(
        gravity                = planet.gravity,
        heat_capacity          = atmosphere.heat_capacity,
        stefan_boltzmann       = atmosphere.stefan_boltzmann,
        solar_constant         = planet.solar_constant,
        dry_air_molar_mass     = atmosphere.mol_mass_dry_air / 1000,
        water_molar_mass       = atmosphere.mol_mass_vapor / 1000,
        dry_air_gas_constant   = atmosphere.R_dry,
        universal_gas_constant = atmosphere.R_gas,
    )
end

@inline function speedy_column_geometry(model)
    geom = model.geometry
    return ColumnGrid(geom.σ_levels_full, geom.σ_levels_half, geom.σ_levels_thick)
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, variables,
                                                                   scheme::AnalyticBandLongwave{NF},
                                                                   model) where NF
    time_stepping = model.time_stepping
    T_all    = SpeedyWeather.get_prognostic_step(variables.grid.temperature, time_stepping, scheme)
    q_all    = SpeedyWeather.get_prognostic_step(variables.grid.humidity, time_stepping, scheme)
    dTdt_all = SpeedyWeather.get_tendency_step(variables.tendencies.grid.temperature, time_stepping, scheme)
    sst_all  = SpeedyWeather.get_prognostic_step(variables.prognostic.ocean.sea_surface_temperature,
                                                 time_stepping, scheme)

    T    = @view T_all[ij, :]
    q    = @view q_all[ij, :]
    Φ    = @view variables.dynamics.geopotential[ij, :]
    temperature_tendency = @view dTdt_all[ij, :]
    pˢ   = variables.parameterizations.surface_pressure[ij]            # [Pa]

    CO₂ = let prog = variables.prognostic
        if hasproperty(prog, :greenhouse_gases) && haskey(prog.greenhouse_gases, :co2)
            NF(prog.greenhouse_gases.co2[])
        else
            NF(DEFAULT_CO₂)
        end
    end

    profile  = AtmosphereProfile(temperature = T, humidity = q,
                                 geopotential = Φ, surface_pressure = pˢ,
                                 CO₂ = CO₂)

    geometry = speedy_column_geometry(model)
    surface  = SurfaceState{NF}(
        sea_surface_temperature  = sst_all[ij],
        land_surface_temperature = variables.prognostic.land.soil_temperature[ij, 1],
        land_fraction            = model.land_sea_mask.land_fraction[ij],
    )
    constants = speedy_physical_constants(model)
    diagnostics = LongwaveDiagnostics{NF}()

    solve_longwave!(temperature_tendency, diagnostics, scheme, profile, geometry, surface, constants)

    variables.parameterizations.outgoing_longwave[ij]         = diagnostics.outgoing_longwave
    variables.parameterizations.surface_longwave_down[ij]     = diagnostics.surface_longwave_down
    variables.parameterizations.surface_longwave_up[ij]       = diagnostics.surface_longwave_up
    variables.parameterizations.ocean.surface_longwave_up[ij] = diagnostics.ocean_surface_longwave_up
    variables.parameterizations.land.surface_longwave_up[ij]  = diagnostics.land_surface_longwave_up

    return nothing
end
