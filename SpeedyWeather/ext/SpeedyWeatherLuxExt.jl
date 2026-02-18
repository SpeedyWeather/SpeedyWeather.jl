module SpeedyWeatherLuxExt

using SpeedyWeather, Lux, NPZ
import Random, Adapt

function SpeedyWeather.LearnedSurfaceRoughness(
        SG::SpectralGrid;
        land_path::String = "z0_land_model_weights.npz",
        kwargs...
    )
    full_land_path = joinpath(@__DIR__, "../input_data/nn_weights", land_path)

    if !isfile(full_land_path)
        @warn "Model weights not found. Ensure paths are correct."
    end

    # Set up Lux NN
    land_nn = Lux.Chain(
        Lux.Dense(13 => 32, Lux.leakyrelu),
        Lux.Dense(32 => 64, Lux.leakyrelu),
        Lux.Dropout(0.2),
        Lux.Dense(64 => 64, Lux.leakyrelu),
        Lux.Dropout(0.1),
        Lux.Dense(64 => 32, Lux.leakyrelu),
        Lux.Dense(32 => 1)
    )
    rng = Random.default_rng()
    land_params, rand_states = Lux.setup(rng, land_nn)
    land_weights = NPZ.npzread(full_land_path)
    load_land_parameters!(land_params, land_weights)

    land_states = Lux.testmode(rand_states)

    return LearnedSurfaceRoughness{SG.NF, typeof(land_nn), typeof(land_params), typeof(land_states)}(;
        land_nn = land_nn,
        land_params = land_params,
        land_states = land_states,
        kwargs...
    )
end

function load_land_parameters!(params, weights::Dict{String, Array{Float32}})
    layer_map = [
        "embed_layer" => :layer_1,
        "layer_1" => :layer_2,
        "layer_2" => :layer_4,
        "layer_3" => :layer_6,
        "output_layer" => :layer_7,
    ]
    for (py_name, lux_sym) in layer_map
        # Check layer names match
        if hasproperty(params, lux_sym)
            lux_layer_params = getproperty(params, lux_sym)
            lux_layer_params.weight .= weights[py_name * ".weight"]
            lux_layer_params.bias .= weights[py_name * ".bias"]
        else
            @warn "Layer $lux_sym not found in Lux params"
        end
    end
    return
end

Adapt.@adapt_structure SpeedyWeather.LearnedSurfaceRoughness
SpeedyWeather.initialize!(::LearnedSurfaceRoughness, ::PrimitiveEquation) = nothing

# Symbolic regression-derived expression for ice-free ocean surface roughness
@inline function ice_free_roughness(x1, x2, x3)
    return (x3 * ((((((-0.096541196 * x3) - 0.30319092) / exp(x3)) - -0.7495353) / exp(x3)) - -0.7495353)) - (-0.040766064 * (((-0.6404533 + x2) * x2) - (x1 - 0.69219255)))
end

@inline function sea_ice_roughness(ℵ)
    return max(1.0f-3, 0.93f-3 * (1.0f0 - ℵ) + (6.05f-3 * exp(-17.0f0 * (ℵ - 0.5f0)^2)))
end

@inline function normalise(a, m, s)
    return (a - m) / s
end

Base.@propagate_inbounds function surface_roughness_land(ij, diagn, progn, scheme::LearnedSurfaceRoughness)
    vₕ = diagn.physics.land.vegetation_high[ij]
    vₗ = diagn.physics.land.vegetation_low[ij]
    vᵦ = 1 - vₕ - vₗ  # bare soil
    g = diagn.grid.geopotential[ij, end]
    sd = progn.land.snow_depth[ij]
    soil_moisture = progn.land.soil_moisture[ij, begin]  # currently top layer
    soil_temperature = progn.land.soil_temperature[ij, end]  # currently bottom layer

    # Modelled interactions
    i₁ = soil_temperature * soil_moisture
    i₂ = soil_temperature * sd
    i₃ = soil_temperature * vₕ
    i₄ = soil_moisture * vₕ
    i₅ = sd * vᵦ
    i₆ = soil_temperature * vᵦ

    # Normalise inputs
    vₕ = normalise(vₕ, scheme.land_input_means[2], scheme.land_input_stds[2])
    vₗ = normalise(vₗ, scheme.land_input_means[3], scheme.land_input_stds[3])
    vᵦ = normalise(vᵦ, scheme.land_input_means[1], scheme.land_input_stds[1])
    g = normalise(g, scheme.land_input_means[4], scheme.land_input_stds[4])
    sd = normalise(sd, scheme.land_input_means[5], scheme.land_input_stds[5])
    soil_moisture = normalise(soil_moisture, scheme.land_input_means[12], scheme.land_input_stds[12])
    soil_temperature = normalise(soil_temperature, scheme.land_input_means[7], scheme.land_input_stds[7])

    # Normalise interaction terms
    i₁ = normalise(i₁, scheme.land_input_means[11], scheme.land_input_stds[11])
    i₂ = normalise(i₂, scheme.land_input_means[10], scheme.land_input_stds[10])
    i₃ = normalise(i₃, scheme.land_input_means[9], scheme.land_input_stds[9])
    i₄ = normalise(i₄, scheme.land_input_means[13], scheme.land_input_stds[13])
    i₅ = normalise(i₅, scheme.land_input_means[6], scheme.land_input_stds[6])
    i₆ = normalise(i₆, scheme.land_input_means[8], scheme.land_input_stds[8])

    scheme.land_input_buffer[:] .= (vᵦ, vₕ, vₗ, g, sd, i₅, soil_temperature, i₆, i₃, i₂, i₁, soil_moisture, i₄)

    prediction, _ = Lux.apply(scheme.land_nn, scheme.land_input_buffer, scheme.land_params, scheme.land_states)
    log_surface_roughness = (prediction[1] * scheme.land_output_std) + scheme.land_output_mean
    surface_roughness = exp(log_surface_roughness)
    return surface_roughness
end

Base.@propagate_inbounds function surface_roughness_ocean(ij, diagn, progn, scheme::LearnedSurfaceRoughness)
    surface = diagn.nlayers
    ℵ = progn.ocean.sea_ice_concentration[ij]
    Uₛ = diagn.grid.u_grid[ij, surface]
    Vₛ = diagn.grid.v_grid[ij, surface]
    UVₛ = diagn.physics.surface_wind_speed[ij]

    Uₛ = normalise(Uₛ, scheme.ocean_input_means[1], scheme.ocean_input_stds[1])
    Vₛ = normalise(Vₛ, scheme.ocean_input_means[2], scheme.ocean_input_stds[2])
    UVₛ = normalise(UVₛ, scheme.ocean_input_means[3], scheme.ocean_input_stds[3])

    log_ocean_roughness = ice_free_roughness(Uₛ, Vₛ, UVₛ) * scheme.ocean_output_std + scheme.ocean_output_mean
    ocean_roughness = exp(log_ocean_roughness)
    ℵ_roughness = sea_ice_roughness(ℵ)  # From IFS documentation, CY49R1

    surface_roughness = ℵ * ℵ_roughness + (1 - ℵ) * ocean_roughness
    return surface_roughness
end

Base.@propagate_inbounds function SpeedyWeather.surface_roughness!(ij, diagn, progn, scheme::LearnedSurfaceRoughness, land_sea_mask)
    land_fraction = land_sea_mask.mask[ij]

    # Compute separate ocean and land surface roughness
    diagn.physics.land.surface_roughness[ij] = land_fraction > 0 ? surface_roughness_land(ij, diagn, progn, scheme) : zero(land_fraction)
    diagn.physics.ocean.surface_roughness[ij] = land_fraction < 1 ? surface_roughness_ocean(ij, diagn, progn, scheme) : zero(land_fraction)

    # Blend the two via arithmetic average
    diagn.physics.surface_roughness[ij] = land_fraction * diagn.physics.land.surface_roughness[ij] + (1 - land_fraction) * diagn.physics.ocean.surface_roughness[ij]
    return nothing
end

end
