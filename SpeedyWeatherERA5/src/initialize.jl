# Stage 3: put the prepared fields into a SpeedyWeather simulation (only existing SpeedyWeather API: set!).

_step(var) = SpeedyWeatherInternals.ArrayDimensions.hastime(var) ? 1 : nothing    # prescribed SST/ice have no step dimension

"""
Set `sim` (a freshly `initialize!`d simulation of a PrimitiveWetModel) from the prepared ERA5 fields `ic` (see `prepare_era5`):
atmosphere, SST and sea-ice concentration. Works with dynamic (`SlabOcean`, `ThermodynamicSeaIce`) and prescribed
(`PrescribedOcean`, `PrescribedSeaIce`) ocean and ice. Land (soil temperature, soil moisture, snow) is not touched."""
function set_era5!(sim, ic)
    model = sim.model
    set!(sim; temperature = ic.temperature, humidity = ic.humidity, u = ic.u, v = ic.v, pressure = ic.pressure)
    ocean = sim.variables.prognostic.ocean
    # SST: cap at the freezing temperature if there is a sea-ice model with one, as SlabOcean's own initialisation does
    sst = hasproperty(model, :sea_ice) && hasproperty(model.sea_ice, :freezing_temperature) ? max.(ic.sst, model.sea_ice.freezing_temperature) : ic.sst
    set!(sim; sea_surface_temperature = sst, namespace = :ocean, step = _step(ocean.sea_surface_temperature))
    set!(sim; sea_ice_concentration = ic.sea_ice, namespace = :ocean, step = _step(ocean.sea_ice_concentration))
    land_temperature = hasproperty(model.ocean, :land_temperature) ? model.ocean.land_temperature : 285
    SpeedyWeather.mask!(ocean.sea_surface_temperature, model.land_sea_mask, :land; masked_value = land_temperature)
    return sim
end

"""
`sim = era5_initialize!(model, dt; cache, kwargs...)`: `initialize!(model; time = dt)`, fetch ERA5 at `dt` if it is not cached,
prepare it for `model` (`kwargs` go to `prepare_era5`) and set the state. `model` must use `Grid = FullGaussianGrid`."""
function era5_initialize!(model, dt::DateTime; cache = default_cache(), kwargs...)
    path = get_era5(dt; cache)
    sim = initialize!(model; time = dt)
    ic = prepare_era5(model, path; kwargs...)
    set_era5!(sim, ic)
    return sim, ic
end
