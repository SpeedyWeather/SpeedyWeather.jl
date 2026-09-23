# NumericalRadiation's ClearSkyEcCKDRadiation (clear-sky ecCKD gas optics, longwave and shortwave
# solved from one gas-optics evaluation) as a single SpeedyWeather `radiation` component:
# PrimitiveWetModel(spectral_grid; radiation = ClearSkyEcCKDRadiation(spectral_grid)). The scheme type
# is NumericalRadiation's own; this file adds the methods SpeedyWeather calls.

"""$(TYPEDSIGNATURES)
An `ClearSkyEcCKDRadiation` in the number format of the spectral grid from tabulated `gas_optics`;
`kwargs` (`ozone`, `mole_fractions`, `surface_emissivity`) as for `NumericalRadiation.ClearSkyEcCKDRadiation`.
CO₂ follows the model's `greenhouse_gases.co2` when present, else the scheme's prescribed
`mole_fractions.co2` (280 ppm by default)."""
ClearSkyEcCKDRadiation(SG::SpeedyWeather.SpectralGrid, gas_optics::EcCKDTabulatedGasOpticsModel; kwargs...) =
    ClearSkyEcCKDRadiation{SG.NF}(gas_optics; kwargs...)

"""$(TYPEDSIGNATURES)
Load the reference ecCKD `model_pair` (`"32x32"` by default, or e.g. `"64x96"`) with the gases
`names`; the tables are read with NCDatasets, a dependency of SpeedyWeather."""
ClearSkyEcCKDRadiation(SG::SpeedyWeather.SpectralGrid, model_pair::Union{AbstractString, Symbol} = "32x32";
               names = (:composite, :h2o, :o3, :co2), kwargs...) =
    ClearSkyEcCKDRadiation(SG.NF, model_pair; names, kwargs...)

SpeedyWeather.initialize!(::ClearSkyEcCKDRadiation, ::SpeedyWeather.PrimitiveEquation) = nothing

number_of_gases(rad::ClearSkyEcCKDRadiation) = length(NumericalRadiation.gas_names(rad.gas_optics))

function SpeedyWeather.variables(rad::ClearSkyEcCKDRadiation{NF}, model::SpeedyWeather.AbstractModel) where NF
    nlayers = SpeedyWeather.get_nlayers(model)
    ng_lw = length(rad.gas_optics.longwave_weights)
    ng_sw = length(rad.gas_optics.shortwave_weights)
    ngas = number_of_gases(rad)
    PV = SpeedyWeather.ParameterizationVariable
    ns = :ecckd
    layers = SpeedyWeather.GridXYZ()
    interfaces = SpeedyWeather.Grid3D(n = nlayers + 1)
    return (
        # standard shortwave and longwave diagnostics (as any AbstractShortwave/AbstractLongwave)
        SpeedyWeather.variables(SpeedyWeather.TransparentShortwave())...,
        SpeedyWeather.variables(SpeedyWeather.UniformCooling{NF}())...,
        # column state
        PV(:pressure_layers, layers, desc = "Layer pressure", units = "Pa", namespace = ns),
        PV(:pressure_interfaces, interfaces, desc = "Interface pressure", units = "Pa", namespace = ns),
        PV(:temperature_interfaces, interfaces, desc = "Interface temperature", units = "K", namespace = ns),
        PV(:gas_amounts, SpeedyWeather.Grid4D(n = ngas), desc = "Layer gas amounts", units = "mol/m^2", namespace = ns),
        # optical properties per g point
        PV(:longwave_optical_depth, SpeedyWeather.Grid4D(n = ng_lw), desc = "Longwave optical depth", units = "1", namespace = ns),
        PV(:longwave_source, SpeedyWeather.Grid4D(n = ng_lw), desc = "Longwave layer Planck source", units = "W/m^2", namespace = ns),
        PV(:longwave_source_top, SpeedyWeather.Grid4D(n = ng_lw), desc = "Longwave Planck source at layer top", units = "W/m^2", namespace = ns),
        PV(:longwave_source_bottom, SpeedyWeather.Grid4D(n = ng_lw), desc = "Longwave Planck source at layer bottom", units = "W/m^2", namespace = ns),
        PV(:shortwave_optical_depth, SpeedyWeather.Grid4D(n = ng_sw), desc = "Shortwave absorption optical depth", units = "1", namespace = ns),
        PV(:shortwave_rayleigh_optical_depth, SpeedyWeather.Grid4D(n = ng_sw), desc = "Shortwave Rayleigh optical depth", units = "1", namespace = ns),
        PV(:shortwave_scattering_asymmetry, SpeedyWeather.Grid4D(n = ng_sw), desc = "Shortwave scattering asymmetry", units = "1", namespace = ns),
        # fluxes and surface emission
        PV(:longwave_up, interfaces, desc = "Upward longwave flux", units = "W/m^2", namespace = ns),
        PV(:longwave_down, interfaces, desc = "Downward longwave flux", units = "W/m^2", namespace = ns),
        PV(:shortwave_up, interfaces, desc = "Upward shortwave flux", units = "W/m^2", namespace = ns),
        PV(:shortwave_down, interfaces, desc = "Downward shortwave flux", units = "W/m^2", namespace = ns),
        PV(:surface_emission, SpeedyWeather.Grid3D(n = ng_lw), desc = "Surface longwave emission per g point", units = "W/m^2", namespace = ns),
        # shortwave adding-method work arrays, the fields of ShortwaveColumnScratch; all with the
        # interface length so that the seven column views share one type (the struct has one
        # array-type parameter), the layer ones use the first nlayers entries
        PV(:sw_reflectance, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_transmittance, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_direct_reflectance, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_direct_diffuse_transmittance, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_direct_flux, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_stack_albedo, interfaces, desc = "Shortwave work array", namespace = ns),
        PV(:sw_source, interfaces, desc = "Shortwave work array", namespace = ns),
    )
end

# column views: (npoints, n) -> vector, (npoints, nlayers, ng) -> (ng, nlayers) matrix
@inline column(field, ij) = view(field, ij, :)
@inline column_gpoints(field, ij) = PermutedDimsArray(view(field, ij, :, :), (2, 1))

"""$(TYPEDSIGNATURES)
Interface temperatures from layer temperatures `T` at pressures `p`: linear in
pressure between layer centres, the layer-1 temperature at the top of the
atmosphere, and the air temperature extrapolated in pressure from the two lowest
layers at the surface. The skin temperature is deliberately *not* used for the
lowest half level: with a lowest layer 100 hPa thick, that would make the whole
layer radiate downward at skin temperature and feed back onto the surface
(observed to run a land surface away to 400 K within a day). The surface itself
emits at its skin temperature through the boundary condition."""
@inline function interface_temperatures!(T_half, T, p, p_half)
    nlayers = length(T)
    T_half[1] = T[1]
    for k in 2:nlayers
        weight = (p_half[k] - p[k - 1]) / (p[k] - p[k - 1])
        T_half[k] = T[k - 1] + weight * (T[k] - T[k - 1])
    end
    if nlayers > 1
        weight = (p_half[nlayers + 1] - p[nlayers]) / (p[nlayers] - p[nlayers - 1])
        T_half[nlayers + 1] = T[nlayers] + weight * (T[nlayers] - T[nlayers - 1])
    else
        T_half[nlayers + 1] = T[nlayers]
    end
    return T_half
end

@inline mole_fraction(x::Number, p) = x
@inline mole_fraction(f, p) = f(p)

# amount [mol/m²] of gas `name` in a layer with `dry` mol/m² of dry air, `h2o` mol/m² water vapour
@inline gas_amount(::Val{:composite}, rad, dry, h2o, co2, p) = dry
@inline gas_amount(::Val{:h2o}, rad, dry, h2o, co2, p) = h2o
@inline gas_amount(::Val{:co2}, rad, dry, h2o, co2, p) = mole_fraction(co2, p) * dry
@inline gas_amount(::Val{name}, rad, dry, h2o, co2, p) where name =
    mole_fraction(getproperty(rad.mole_fractions, name), p) * dry

"""$(TYPEDSIGNATURES)
Fill `amounts` (shape `(nlayers, ngas)`, gas order `names`) with layer molar
amounts [mol/m²] from specific humidity `q`, layer pressures `p`, interface
pressures `p_half`, the CO₂ mole fraction `co2` (a number, or a function of pressure) and the host's `constants`
(gravity and the molar masses of dry air and water)."""
@generated function gas_amounts!(amounts, ::Val{names}, rad, q, p, p_half, co2, constants) where names
    assignments = [:(amounts[k, $j] = gas_amount(Val($(QuoteNode(name))), rad, dry, h2o, co2, p[k]))
                   for (j, name) in enumerate(names)]
    return quote
        Base.@_propagate_inbounds_meta
        (; gravity, dry_air_molar_mass, water_molar_mass) = constants
        for k in eachindex(q)
            moist_mass = (p_half[k + 1] - p_half[k]) / gravity     # kg/m² of moist air
            water_mass = q[k] * moist_mass
            dry = (moist_mass - water_mass) / dry_air_molar_mass
            h2o = water_mass / water_molar_mass
            $(assignments...)
        end
        return amounts
    end
end

@inline gas_views(amounts, ::Val{names}) where names =
    NamedTuple{names}(ntuple(j -> view(amounts, :, j), Val(length(names))))

# The column update is split into its stages; every stage works on views into the
# `vars.parameterizations.ecckd` work arrays (W) for column ij and allocates nothing.
Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars,
                                                                   rad::ClearSkyEcCKDRadiation{NF},
                                                                   model) where NF
    W = vars.parameterizations.ecckd
    time_stepping = model.time_stepping
    T_all = SpeedyWeather.get_prognostic_step(vars.grid.temperature, time_stepping, rad)
    q_all = SpeedyWeather.get_prognostic_step(vars.grid.humidity, time_stepping, rad)
    dTdt  = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, time_stepping, rad)
    T = @view T_all[ij, :]
    q = @view q_all[ij, :]
    pₛ = vars.parameterizations.surface_pressure[ij]

    surface = ecckd_surface_state(ij, vars, rad, model)
    atmosphere = ecckd_column_atmosphere!(ij, W, rad, T, q, pₛ, surface, vars, model)
    longwave, shortwave = ecckd_column_optics(ij, W, rad)
    optical_properties!(longwave, shortwave, rad.gas_optics, atmosphere)

    fluxes = RadiativeFluxes(longwave_up = column(W.longwave_up, ij),
                             longwave_down = column(W.longwave_down, ij),
                             shortwave_up = column(W.shortwave_up, ij),
                             shortwave_down = column(W.shortwave_down, ij))
    ecckd_longwave!(ij, vars, fluxes, longwave, atmosphere, rad, surface)
    ecckd_shortwave!(ij, vars, fluxes, shortwave, atmosphere, rad, surface, model)
    ecckd_heating!(ij, dTdt, fluxes, pₛ, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Lower-boundary state of column `ij`: land fraction, sea-surface and land-surface
temperature and their blend, ocean and land albedo and their blend, cosine of the
solar zenith angle."""
Base.@propagate_inbounds function ecckd_surface_state(ij, vars, rad::ClearSkyEcCKDRadiation, model)
    sst_all = SpeedyWeather.get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature,
                                                model.time_stepping, rad)
    land_fraction = model.land_sea_mask.land_fraction[ij]
    sea_surface_temperature = sst_all[ij]
    land_surface_temperature = vars.prognostic.land.soil_temperature[ij, 1]
    albedo_ocean = vars.parameterizations.ocean.albedo[ij]
    albedo_land = vars.parameterizations.land.albedo[ij]
    return (;
        land_fraction,
        sea_surface_temperature,
        land_surface_temperature,
        temperature = (1 - land_fraction) * sea_surface_temperature + land_fraction * land_surface_temperature,
        albedo_ocean,
        albedo_land,
        albedo = (1 - land_fraction) * albedo_ocean + land_fraction * albedo_land,
        cos_zenith = vars.parameterizations.cos_zenith[ij],
    )
end

"""$(TYPEDSIGNATURES)
Fill the column's layer and interface pressures, interface temperatures and gas
amounts in the work arrays `W` and return them as a `ColumnAtmosphere`. The
blended skin temperature is carried in `surface.temperature` for the solvers'
boundary conditions only. CO₂ is
taken from the model's greenhouse gases when present, else from the scheme's prescribed
`mole_fractions.co2`."""
Base.@propagate_inbounds function ecckd_column_atmosphere!(ij, W, rad::ClearSkyEcCKDRadiation{NF}, T, q, pₛ,
                                                           surface, vars, model) where NF
    nlayers = length(T)
    coordinates = model.geometry.vertical_coordinates

    p = column(W.pressure_layers, ij)
    p_half = column(W.pressure_interfaces, ij)
    T_half = column(W.temperature_interfaces, ij)
    for k in 1:nlayers
        p[k] = SpeedyWeather.pressure(k, pₛ, coordinates)
        p_half[k] = SpeedyWeather.pressure_half(k, pₛ, coordinates)
    end
    p_half[nlayers + 1] = SpeedyWeather.pressure_half(nlayers + 1, pₛ, coordinates)
    interface_temperatures!(T_half, T, p, p_half)

    co2 = let prognostic = vars.prognostic
        if hasproperty(prognostic, :greenhouse_gases) && haskey(prognostic.greenhouse_gases, :co2)
            NF(prognostic.greenhouse_gases.co2[]) * NF(1e-6)        # [ppm] → mole fraction
        else
            rad.mole_fractions.co2                                  # number or function of pressure
        end
    end
    names = Val(NumericalRadiation.gas_names(rad.gas_optics))
    amounts = @view W.gas_amounts[ij, :, :]
    constants = speedy_physical_constants(model)
    gas_amounts!(amounts, names, rad, q, p, p_half, co2, constants)

    return ColumnAtmosphere(
        pressure_layers = p, pressure_interfaces = p_half,
        temperature_layers = T, temperature_interfaces = T_half,
        gases = gas_views(amounts, names),
        surface = (; temperature = surface.temperature),
        geometry = (; cos_zenith = surface.cos_zenith),
        constants = constants,
    )
end

"""$(TYPEDSIGNATURES)
`LongwaveOptics` and `ShortwaveOptics` of column `ij` as `(ng, nlayers)` views
into the work arrays `W`, with the gas-optics model's spectral weights."""
Base.@propagate_inbounds function ecckd_column_optics(ij, W, rad::ClearSkyEcCKDRadiation)
    (; gas_optics) = rad
    longwave = LongwaveOptics(column_gpoints(W.longwave_optical_depth, ij),
                              column_gpoints(W.longwave_source, ij);
                              source_top = column_gpoints(W.longwave_source_top, ij),
                              source_bottom = column_gpoints(W.longwave_source_bottom, ij),
                              weights = gas_optics.longwave_weights)
    shortwave = ShortwaveOptics(column_gpoints(W.shortwave_optical_depth, ij);
                                rayleigh_optical_depth = column_gpoints(W.shortwave_rayleigh_optical_depth, ij),
                                scattering_asymmetry = column_gpoints(W.shortwave_scattering_asymmetry, ij),
                                weights = gas_optics.shortwave_weights)
    return longwave, shortwave
end

"""$(TYPEDSIGNATURES)
Longwave stream: spectral surface emission over ocean and land blended by land
fraction, clear-sky transfer, and the longwave diagnostics of column `ij`."""
Base.@propagate_inbounds function ecckd_longwave!(ij, vars, fluxes, longwave, atmosphere,
                                                  rad::ClearSkyEcCKDRadiation{NF}, surface) where NF
    (; gas_optics) = rad
    W = vars.parameterizations.ecckd
    nlayers = length(atmosphere.temperature_layers)
    f = surface.land_fraction

    # per-g-point surface emission over ocean and land, evaluated lazily from the
    # source table and blended by land fraction into the column's work vector
    emission = column(W.surface_emission, ij)
    ocean = TabulatedSurfaceEmission(gas_optics, surface.sea_surface_temperature;
                                     emissivity = rad.surface_emissivity)
    land = TabulatedSurfaceEmission(gas_optics, surface.land_surface_temperature;
                                    emissivity = rad.surface_emissivity)
    weights = gas_optics.longwave_weights
    up_ocean = zero(NF)             # broadband emission over ocean and land for their models
    up_land = zero(NF)
    for ig in eachindex(weights)
        eₒ, eₗ = ocean[ig], land[ig]
        up_ocean += weights[ig] * eₒ
        up_land += weights[ig] * eₗ
        emission[ig] = (1 - f) * eₒ + f * eₗ
    end

    radiative_fluxes!(fluxes, CloudlessLongwave(), longwave, atmosphere,
                      LongwaveBoundaryConditions(surface_longwave_up = emission))

    vars.parameterizations.outgoing_longwave[ij] = fluxes.longwave_up[1]
    vars.parameterizations.surface_longwave_down[ij] = fluxes.longwave_down[nlayers + 1]
    vars.parameterizations.surface_longwave_up[ij] = fluxes.longwave_up[nlayers + 1]
    vars.parameterizations.ocean.surface_longwave_up[ij] = up_ocean
    vars.parameterizations.land.surface_longwave_up[ij] = up_land
    return nothing
end

"""$(TYPEDSIGNATURES)
Shortwave stream of column `ij`: clear-sky transfer with Rayleigh scattering for
the sunlit hemisphere (zero fluxes at night) and the shortwave diagnostics."""
Base.@propagate_inbounds function ecckd_shortwave!(ij, vars, fluxes, shortwave, atmosphere,
                                                   rad::ClearSkyEcCKDRadiation{NF}, surface, model) where NF
    W = vars.parameterizations.ecckd
    nlayers = length(atmosphere.temperature_layers)
    (; cos_zenith, albedo, albedo_ocean, albedo_land) = surface

    if cos_zenith > 0
        scratch = ShortwaveColumnScratch(view(W.sw_reflectance, ij, 1:nlayers),
                                         view(W.sw_transmittance, ij, 1:nlayers),
                                         view(W.sw_direct_reflectance, ij, 1:nlayers),
                                         view(W.sw_direct_diffuse_transmittance, ij, 1:nlayers),
                                         view(W.sw_direct_flux, ij, 1:nlayers),
                                         view(W.sw_stack_albedo, ij, 1:(nlayers + 1)),
                                         view(W.sw_source, ij, 1:(nlayers + 1)))
        boundary = ShortwaveBoundaryConditions(
            toa_shortwave_down = model.planet.solar_constant * cos_zenith,
            surface_albedo = albedo)
        radiative_fluxes!(fluxes, CloudlessShortwave(), shortwave, atmosphere, boundary, scratch)
    else
        for k in 1:(nlayers + 1)
            fluxes.shortwave_up[k] = zero(NF)
            fluxes.shortwave_down[k] = zero(NF)
        end
    end

    surface_down = fluxes.shortwave_down[nlayers + 1]
    vars.parameterizations.outgoing_shortwave[ij] = fluxes.shortwave_up[1]
    vars.parameterizations.surface_shortwave_down[ij] = surface_down
    vars.parameterizations.ocean.surface_shortwave_down[ij] = surface_down
    vars.parameterizations.land.surface_shortwave_down[ij] = surface_down
    vars.parameterizations.surface_shortwave_up[ij] = fluxes.shortwave_up[nlayers + 1]
    vars.parameterizations.ocean.surface_shortwave_up[ij] = albedo_ocean * surface_down
    vars.parameterizations.land.surface_shortwave_up[ij] = albedo_land * surface_down
    vars.parameterizations.albedo[ij] = albedo
    return nothing
end

"""$(TYPEDSIGNATURES)
Add the net (longwave + shortwave) flux convergence of every layer of column `ij`
to the temperature tendency `dTdt`, in K/s, via SpeedyWeather's `flux_to_tendency`."""
Base.@propagate_inbounds function ecckd_heating!(ij, dTdt, fluxes, pₛ, model)
    cₚ = model.atmosphere.heat_capacity
    nlayers = length(fluxes.longwave_up) - 1
    net(k) = fluxes.longwave_down[k] - fluxes.longwave_up[k] +
             fluxes.shortwave_down[k] - fluxes.shortwave_up[k]
    for k in 1:nlayers
        dTdt[ij, k] += SpeedyWeather.flux_to_tendency((net(k) - net(k + 1)) / cₚ, pₛ, k, model)
    end
    return nothing
end
