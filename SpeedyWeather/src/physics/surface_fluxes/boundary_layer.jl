abstract type AbstractBoundaryLayer <: AbstractParameterization end

variables(::AbstractBoundaryLayer) = (
    DiagnosticVariable(name = :boundary_layer_drag, dims = Grid2D(), desc = "Boundary layer drag coefficient", units = "1"),
)

export ConstantDrag
"""Constant boundary layer drag coefficient. Fields are $(TYPEDFIELDS)"""
@kwdef struct ConstantDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] Constant drag coefficient [1]"
    drag::NF = 1.0e-3
end

Adapt.@adapt_structure ConstantDrag
ConstantDrag(SG::SpectralGrid; kwargs...) = ConstantDrag{SG.NF}(; kwargs...)
initialize!(::ConstantDrag, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, diagn, progn, drag::ConstantDrag, model) =
    boundary_layer_drag!(ij, diagn, drag)

@propagate_inbounds function boundary_layer_drag!(ij, diagn, scheme::ConstantDrag)
    return diagn.physics.boundary_layer_drag[ij] = scheme.drag
end

abstract type AbstractSurfaceRoughness <: AbstractParameterization end

export ConstantSurfaceRoughness
@kwdef struct ConstantSurfaceRoughness{NF} <: AbstractSurfaceRoughness
    "[OPTION] constant roughness length over land [m]"
    roughness_length_land::NF = 0.5

    "[OPTION] constant roughness length over ocean [m]"
    roughness_length_ocean::NF = 1.0e-4
end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid, kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing

export LearnedSurfaceRoughness
@kwdef struct LearnedSurfaceRoughness{NF, LM, LNN, LP, LS} <: AbstractSurfaceRoughness
    # Ocean normalisation parameters
    ocean_input_means::Vector{NF} = Float32[0.19490805, 0.11980359, 7.7569385]
    ocean_input_stds::Vector{NF} = Float32[6.6221075, 5.4018555, 3.5962458]
    ocean_output_mean::NF = -9.1571865
    ocean_output_std::NF = 1.0090618

    # Land normalisation parameters
    land_input_means::Vector{NF} = Float32[
        7.4566591e-1, 1.0025085e-1, 1.5397815e-1, 1.6788273e+4,
        6.3441253e+0, 6.3426518e+0, 2.5690454e+2, 1.8242287e+2,
        2.9126474e+1, 1.4939866e+3, 5.4704014e+1, 2.1826939e-1,
        3.3305537e-2,
    ]
    land_input_stds::Vector{NF} = Float32[
        4.23785776e-1, 2.76675612e-1, 3.28937173e-1, 1.16441348e+4,
        4.8089776e+0, 4.81200314e+0, 3.01283646e+1, 1.04803604e+2,
        8.05910187e+1, 1.13503149e+3, 3.34448814e+1, 1.26479045e-1,
        9.95290726e-2,
    ]
    land_output_mean::NF = -5.031811
    land_output_std::NF = 2.4447718

    land_input_buffer::Vector{Float32} = zeros(Float32, 13)

    land_model::LM
    land_nn::LNN
    land_params::LP
    land_states::LS
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

function LearnedSurfaceRoughness(
        SG::SpectralGrid;
        land_path::String = "z0_land_model_weights.npz",
        kwargs...
    )
    full_land_path = joinpath(@__DIR__, "../../../input_data/nn_weights", land_path)

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
    land_model = (u, p, s) -> first(Lux.apply(land_nn, u, p, s))

    return LearnedSurfaceRoughness{SG.NF, typeof(land_model), typeof(land_nn), typeof(land_params), typeof(land_states)}(;
        land_model = land_model,
        land_nn = land_nn,
        land_params = land_params,
        land_states = land_states,
        kwargs...
    )
end

Adapt.@adapt_structure LearnedSurfaceRoughness
initialize!(::LearnedSurfaceRoughness, ::PrimitiveEquation) = nothing
@propagate_inbounds function surface_roughness_land(ij, scheme::ConstantSurfaceRoughness, diagn, progn)
    return scheme.roughness_length_land
end

@propagate_inbounds function surface_roughness_ocean(ij, scheme::ConstantSurfaceRoughness, diagn, progn)
    return scheme.roughness_length_ocean
end

@propagate_inbounds function surface_roughness_land(ij, scheme::LearnedSurfaceRoughness, diagn, progn)
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
    normalise(a, m, s) = (a - m) / s
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

    prediction = scheme.land_model(scheme.land_input_buffer, scheme.land_params, scheme.land_states)
    log_surface_roughness = (prediction[1] * scheme.land_output_std) + scheme.land_output_mean
    surface_roughness = exp(log_surface_roughness)
    return surface_roughness
end

@propagate_inbounds function surface_roughness_ocean(ij, scheme::LearnedSurfaceRoughness, diagn, progn)
    surface = diagn.nlayers
    ℵ = progn.ocean.sea_ice_concentration[ij]
    Uₛ = diagn.grid.u_grid[ij, surface]
    Vₛ = diagn.grid.v_grid[ij, surface]
    UVₛ = diagn.physics.surface_wind_speed[ij]

    normalise(a, m, s) = (a - m) / s
    Uₛ = normalise(Uₛ, scheme.ocean_input_means[1], scheme.ocean_input_stds[1])
    Vₛ = normalise(Vₛ, scheme.ocean_input_means[2], scheme.ocean_input_stds[2])
    UVₛ = normalise(UVₛ, scheme.ocean_input_means[3], scheme.ocean_input_stds[3])

    # Symbolic regression-derived expression for ice-free ocean surface roughness
    ice_free_roughness(x1, x2, x3) = max(abs(-1.2154098f0 - (((0.22197704f0 * min(x2, x3)) + max(exp(0.7497865f0 - x3), x1)) * -0.19615653f0)) - 0.73932385f0, -0.40111288f0) + (x3 / 1.2592316f0)
    log_ocean_roughness = ice_free_roughness(Uₛ, Vₛ, UVₛ) * scheme.ocean_output_std + scheme.ocean_output_mean
    ocean_roughness = exp(log_ocean_roughness)
    sea_ice_roughness = max(1.0f-3, 0.93f-3 * (1.0f0 - ℵ) + (6.05f-3 * exp(-17.0f0 * (ℵ - 0.5f0)^2)))  # From IFS documentation, CY49R1

    surface_roughness = ℵ * sea_ice_roughness + (1 - ℵ) * ocean_roughness
    return surface_roughness
end

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@kwdef struct BulkRichardsonDrag{NF, SR} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    von_Karman::NF = 0.4

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    critical_Richardson::NF = 10

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    drag_min::NF = 1.0e-5

    surface_roughness::SR
end

Adapt.@adapt_structure BulkRichardsonDrag
function BulkRichardsonDrag(
        SG::SpectralGrid;
        surface_roughness = ConstantSurfaceRoughness(SG),
        kwargs...
    )
    SR = typeof(surface_roughness)
    return BulkRichardsonDrag{SG.NF, SR}(; surface_roughness, kwargs...)
end
function initialize!(drag::BulkRichardsonDrag, model::PrimitiveEquation)
    initialize!(drag.surface_roughness, model)
    return nothing
end

# function barrier
@propagate_inbounds function parameterization!(ij, diagn, progn, drag::BulkRichardsonDrag, model)
    boundary_layer_drag!(ij, diagn, progn, drag, model.land_sea_mask, model.atmosphere, model.planet, model.orography)
    return boundary_layer_drag!(ij, diagn, progn, drag, model.land_sea_mask, model.atmosphere, model.planet, model.orography)
end

@propagate_inbounds function boundary_layer_drag!(
        ij,
        diagn,
        progn,
        drag::BulkRichardsonDrag,
        land_sea_mask,
        atmosphere,
        planet,
        orography,
    )
    land_fraction = land_sea_mask.mask[ij]

    # Height z [m] of lowermost layer above ground
    surface = diagn.nlayers
    (; gravity) = planet
    z = diagn.grid.geopotential[ij, surface] / gravity - orography.orography[ij]

    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    κ = drag.von_Karman

    # Calculate land and ocean roughness lengths
    z₀_land = surface_roughness_land(ij, drag.surface_roughness, diagn, progn)
    z₀_ocean = surface_roughness_ocean(ij, drag.surface_roughness, diagn, progn)

    z₀ = land_fraction * z₀_land + (1 - land_fraction) * z₀_ocean

    # should be z > z₀, z=z₀ means an infinitely high drag
    # 0 < z < z₀ doesn't make sense so cap here
    z = max(z, z₀)
    drag_max = (κ / log(z / z₀))^2

    # bulk Richardson number at lowermost layer from Frierson, 2006, eq. (15)
    # they call it Ri_a = Ri here
    ΔΦ₀ = gravity * z     # geopotential high relative to surface
    Ri = bulk_richardson_surface(ij, ΔΦ₀, diagn, atmosphere)
    Ri_c = drag.critical_Richardson
    (; drag_min) = drag

    # clamp to get the cases, eq (12-14)
    # if Ri > Ri_c then C = 0
    # if Ri_c > Ri > 0 then = κ^2/log(z/z₀)^2 * (1-Ri/Ri_c)^2
    # if Ri_c < 0 then κ^2/log(z/z₀)^2
    Ri = clamp(Ri, 0, Ri_c)
    diagn.physics.boundary_layer_drag[ij] = max(drag_min, drag_max * (1 - Ri / Ri_c)^2)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk Richardson number following Frierson, 2006.
For vertical stability in the boundary layer."""
@propagate_inbounds function bulk_richardson_surface(ij, ΔΦ₀, diagn, atmosphere)
    cₚ = atmosphere.heat_capacity
    surface = diagn.nlayers     # surface index = nlayers

    Vₛ = diagn.physics.surface_wind_speed[ij]
    T = diagn.grid.temp_grid_prev[ij, surface]
    q = diagn.grid.humid_grid_prev[ij, surface]
    Tᵥ = virtual_temperature(T, q, atmosphere)

    # bulk Richardson number at lowermost layer N from Frierson, 2006, eq. (15)
    Θ₀ = cₚ * Tᵥ          # virtual dry static energy at surface (z=0)
    Θ₁ = Θ₀ + ΔΦ₀       # virtual dry static energy at first model level (z=z)
    bulk_richardson = ΔΦ₀ * (Θ₁ - Θ₀) / (Θ₀ * Vₛ^2)
    return bulk_richardson
end
