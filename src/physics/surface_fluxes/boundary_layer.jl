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
@kwdef struct LearnedSurfaceRoughness{NF, OM, LM} <: AbstractSurfaceRoughness
    "ONNX Inference Session for Ocean"
    ocean_model::OM

    # Ocean normalisation parameters
    ocean_input_means::Vector{NF} = Float32[0.14675693, 0.24450141, 0.17968568, 7.6465526]
    ocean_input_stds::Vector{NF} = Float32[0.3357911, 6.3831706, 5.53958, 3.6144474]
    ocean_output_mean::NF = -8.76918
    ocean_output_std::NF = 1.2418048

    "ONNX Inference Session for Land"
    land_model::LM

    # Land normalisation parameters
    land_input_means::Vector{NF} = Float32[0.100255094, 0.154690117, 16791.6582, 6.33934355]
    land_input_stds::Vector{NF} = Float32[0.27687377, 0.32951584, 11649.544, 4.8004246]
    land_output_mean::NF = -5.0261893
    land_output_std::NF = 2.4474483
end

function LearnedSurfaceRoughness(
        SG::SpectralGrid;
        ocean_path = "ocean_model/z0_ocean_model_smol.onnx",
        land_path = "land_model/xgb_land",
        kwargs...
    )

    full_ocean_path = joinpath(@__DIR__, "../../../input_data", ocean_path)
    full_land_path = joinpath(@__DIR__, "../../../input_data", land_path)

    if !isfile(full_ocean_path) || !isfile(full_land_path)
        @warn "ONNX models not found. Ensure paths are correct."
    end

    ocean_model = ONNXRunTime.load_inference(full_ocean_path)
    # land_model = ONNXRunTime.load_inference(full_land_path)
    land_model = XGBoost.load(XGBoost.Booster, full_land_path)

    return LearnedSurfaceRoughness{SG.NF, typeof(ocean_model), typeof(land_model)}(;
        ocean_model = ocean_model,
        land_model = land_model,
        kwargs...
    )
end

Adapt.@adapt_structure LearnedSurfaceRoughness
LearnedSurfaceRoughness(SG::SpectralGrid, kwargs...) = LearnedSurfaceRoughness{SG.NF, OM, LM}(; kwargs...)
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
    g = diagn.grid.geopotential[ij, 8]  # surface layer
    sd = progn.land.snow_depth[ij]

    raw_inputs = Float32[vₕ vₗ g sd]
    inputs_norm = Float32.(real.((raw_inputs .- scheme.land_input_means) ./ scheme.land_input_stds))
    # surface_roughness = XGBoost.predict(scheme.land_model, raw_inputs)[1]

    # Symbolic regression replacement expression
    surface_roughness = min(max(0.0045800665, min(inputs_norm[1] / max(inputs_norm[2], 1.4739709), (1.2306644 + inputs_norm[4]) / -0.0021820515)), max(1.2238768 - inputs_norm[2], 0.6677175))
    return Float32(max(surface_roughness, 1.2999999e-3))
end

@propagate_inbounds function surface_roughness_ocean(ij, scheme::LearnedSurfaceRoughness, diagn, progn)
    surface = diagn.nlayers

    # Extract features
    ℵ = progn.ocean.sea_ice_concentration[ij]
    Uₛ = diagn.grid.u_grid[ij, surface]
    Vₛ = diagn.grid.v_grid[ij, surface]
    UVₛ = diagn.physics.surface_wind_speed[ij]

    raw_inputs = [ℵ, Uₛ, Vₛ, UVₛ]
    inputs_norm = Float32.((raw_inputs .- scheme.ocean_input_means) ./ scheme.ocean_input_stds)

    predictors = reshape(inputs_norm, 1, 4)
    input = Dict("x_input" => predictors)

    # Run Inference using the model stored in the struct
    prediction = scheme.ocean_model(input)["linear_5"][1]

    # Post-process
    log_surface_roughness = (prediction * scheme.ocean_output_std) + scheme.ocean_output_mean
    surface_roughness = exp(log_surface_roughness)
    return Float32(max(surface_roughness, 2.4361885e-5))
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

    surface_roughness::SR = ConstantSurfaceRoughness{NF}()
end

Adapt.@adapt_structure BulkRichardsonDrag
function BulkRichardsonDrag(
        SG::SpectralGrid;
        surface_roughness = ConstantSurfaceRoughness(SG), # Default if not provided
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
Calculate the bulk richardson number following Frierson, 2006.
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
