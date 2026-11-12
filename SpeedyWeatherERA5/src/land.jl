# Land: ERA5 soil temperature, soil moisture and snow -> SpeedyWeather's 2-layer bucket land (LandBucketTemperature,
# LandBucketMoisture, SnowModel). The mapping choices are documented in NOTES.md ("Land").

const ERA5_SOIL_BOTTOMS = [0.07, 0.28, 1.00, 2.89]              # m, bottoms of the four ERA5 soil layers (IFS)

"""Row k holds the weights (summing to 1) of the four ERA5 soil layers in model layer k, from the overlap of the depth ranges."""
function soil_layer_weights(thickness)
    model_bottoms = cumsum(thickness)
    model_tops = [0.0; model_bottoms[1:(end - 1)]]
    era_tops = [0.0; ERA5_SOIL_BOTTOMS[1:(end - 1)]]
    W = [max(0.0, min(model_bottoms[k], ERA5_SOIL_BOTTOMS[l]) - max(model_tops[k], era_tops[l])) for k in eachindex(thickness), l in 1:4]
    return W ./ sum(W, dims = 2)
end

"""
Prepare ERA5 land for `model` (a PrimitiveWetModel with `LandModel`): `atmosphere_path` is the file from `get_era5` (land-sea mask,
snow depth), `land_path` the one from `get_era5_land`. Returns `(; soil_temperature, soil_moisture, snow_depth, info)`:
  soil_temperature [K], 2 layers: depth-weighted (layer overlap) means of ERA5's 4 layers on the model's layer thicknesses
  soil_moisture_volumetric [m³/m³], 2 layers: depth-weighted means of ERA5's volumetric soil water, clamped to [0, 1]
  soil_moisture_relative [1], 2 layers: the same divided by the model's `field_capacity`, clamped to [0, 1]
                    (the convention `LandBucketMoisture`'s dynamics document: 1 = field capacity, anything above is removed)
  NOTE SpeedyWeather's own soil-moisture climatology is in the VOLUMETRIC convention (land mean 0.24 against ERA5's 0.22 m³/m³)
  although the dynamics use the relative one, so both are offered, see `set_era5_land!` and NOTES.md.
  snow_depth [m of liquid water equivalent], capped at the model's `snow_depth_cap`
All averaged onto the model grid over ERA5 LAND cells only (lsm > 0.5), cells without land are filled from neighbours
(the caller masks the ocean points as SpeedyWeather does)."""
function prepare_era5_land(model, atmosphere_path, land_path)
    thickness = Float64.(Array(model.land.geometry.layer_thickness))
    @assert length(thickness) == 2 "land layers other than 2 are not supported"
    W = soil_layer_weights(thickness)
    γ = Float64(model.land.thermodynamics.field_capacity)
    cap = Float64(model.land.snow.snow_depth_cap)

    dsa = NCDataset(atmosphere_path); dsl = NCDataset(land_path)
    rd(ds, v) = Float64.(ds[v][:, 2:720])                       # drop pole rows, as for the atmosphere
    land = rd(dsa, "land_sea_mask") .> 0.5
    snow = rd(dsa, "snow_depth")
    stl = [rd(dsl, "soil_temperature_level_$l") for l in 1:4]
    swvl = [rd(dsl, "volumetric_soil_water_layer_$l") for l in 1:4]
    close(dsa); close(dsl)

    landonly(a) = ifelse.(land, a, NaN)
    layer(f, k) = sum(W[k, l] .* f[l] for l in 1:4)             # depth-weighted mean of ERA5 layers for model layer k
    to_model(a) = Matrix(masked_to_model(model, landonly(a)))

    T = cat((to_model(layer(stl, k)) for k in 1:2)...; dims = 3)
    Mvol = cat((to_model(layer(swvl, k)) for k in 1:2)...; dims = 3)
    S = min.(to_model(snow), cap)
    info = (; soil_layer_weights = W, field_capacity = γ, snow_cap = cap,
              relative_moisture_unclamped_max = maximum(Mvol) / γ)
    field(a) = FullGaussianField(Float32.(a), input_as = Matrix)
    return (; soil_temperature = field(T), soil_moisture_volumetric = field(clamp.(Mvol, 0, 1)), soil_moisture_relative = field(clamp.(Mvol ./ γ, 0, 1)),
            snow_depth = field(S), info)
end

# the climatology initialisation leaves both leapfrog steps equal; do the same
function _set_all_steps!(var)
    d = var.data
    n = size(d, ndims(d))
    for s in 2:n
        selectdim(d, ndims(d), s) .= selectdim(d, ndims(d), 1)
    end
    return var
end

"""
Overwrite SpeedyWeather's land climatology (soil temperature, soil moisture, snow = 0) in `sim` by the ERA5 values `land`
from `prepare_era5_land`. `temperature` and `snow` are `true`/`false`; `moisture` is `:relative` (ERA5 / field capacity, the dynamics'
convention), `:volumetric` (raw m³/m³, the convention of SpeedyWeather's own climatology) or `false`. Anything switched off stays
on climatology.
Points that are entirely ocean get SpeedyWeather's own ocean values, as its own initialisation does."""
function set_era5_land!(sim, land; temperature = true, moisture = :relative, snow = true)
    model, lv = sim.model, sim.variables.prognostic.land
    if temperature
        set!(sim; soil_temperature = land.soil_temperature, namespace = :land, step = _step(lv.soil_temperature))
        _set_all_steps!(lv.soil_temperature)
        SpeedyWeather.mask!(lv.soil_temperature, model.land_sea_mask, :ocean; masked_value = model.land.temperature.ocean_temperature)
    end
    if moisture !== false
        @assert moisture in (:relative, :volumetric) "moisture must be :relative, :volumetric or false"
        set!(sim; soil_moisture = moisture == :relative ? land.soil_moisture_relative : land.soil_moisture_volumetric, namespace = :land, step = _step(lv.soil_moisture))
        _set_all_steps!(lv.soil_moisture)
        SpeedyWeather.mask!(lv.soil_moisture, model.land_sea_mask, :ocean; masked_value = model.land.soil_moisture.ocean_moisture)
    end
    if snow
        set!(sim; snow_depth = land.snow_depth, namespace = :land, step = _step(lv.snow_depth))
        _set_all_steps!(lv.snow_depth)
        SpeedyWeather.mask!(lv.snow_depth, model.land_sea_mask, :ocean; masked_value = 0)
    end
    return sim
end
