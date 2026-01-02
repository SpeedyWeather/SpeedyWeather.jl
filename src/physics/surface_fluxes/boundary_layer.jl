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

export BulkRichardsonDrag

"""Boundary layer drag coefficient from the bulk Richardson number,
following Frierson, 2006, Journal of the Atmospheric Sciences.
$(TYPEDFIELDS)"""
@kwdef struct BulkRichardsonDrag{NF} <: AbstractBoundaryLayer
    "[OPTION] von Kármán constant [1]"
    von_Karman::NF = 0.4

    "[OPTION] roughness length over land [m]"
    roughness_length_land::NF = 0.5

    "[OPTION] roughness length over ocean [m]"
    roughness_length_ocean::NF = 1.0e-4

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    critical_Richardson::NF = 10

    "[OPTION] Drag minimum to avoid zero surface fluxes in stable conditions [1]"
    drag_min::NF = 1.0e-5
end

Adapt.@adapt_structure BulkRichardsonDrag
BulkRichardsonDrag(SG::SpectralGrid, kwargs...) = BulkRichardsonDrag{SG.NF}(; kwargs...)
initialize!(::BulkRichardsonDrag, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds parameterization!(ij, diagn, progn, drag::BulkRichardsonDrag, model) =
    boundary_layer_drag!(ij, diagn, progn, drag, model.land_sea_mask, model.atmosphere, model.planet, model.orography)

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
    # z₀ = land_fraction * drag.roughness_length_land + (1 - land_fraction) * drag.roughness_length_ocean
    z₀ = land_fraction * surface_roughness_land(ij, diagn, progn) + (1 - land_fraction) * surface_roughness_ocean(ij, diagn, progn)

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

const OCEAN_ROUGHNESS_MODEL_REF = Ref{ONNXRunTime.InferenceSession}()
const LAND_ROUGHNESS_MODEL_REF = Ref{ONNXRunTime.InferenceSession}()

# 2. A private helper function to avoid repeating code
function _load_model_safe!(ref_container, filename)
    path = joinpath(@__DIR__, "../../../input_data", filename)

    if !isfile(path)
        @warn "ONNX Model not found at $path"
    else
        # Load directly into the container passed as argument
        ref_container[] = ONNXRunTime.load_inference(path)
    end
end

# 3. The SINGLE __init__ function
function __init__()
    # Load Ocean Model
    _load_model_safe!(OCEAN_ROUGHNESS_MODEL_REF, "ocean_model/model.onnx")

    # Load Land Model
    _load_model_safe!(LAND_ROUGHNESS_MODEL_REF, "land_model/rf_z0_land_model.onnx")
end


const O_INPUT_MEANS = Float32[0.14675693, 0.24450141, 0.17968568, 7.6465526]
const O_INPUT_STDS  = Float32[0.3357911, 6.3831706, 5.53958, 3.6144474]
const O_OUTPUT_MEAN = Float32(-8.76918)
const O_OUTPUT_STD  = Float32(1.2418048)

const L_INPUT_MEANS = Float32[0.100255094, 0.154690117, 16791.6582, 6.33934355]
const L_INPUT_STDS  = Float32[0.27687377, 0.32951584, 11649.544, 4.8004246]
const L_OUTPUT_MEAN = Float32(-5.0261893)
const L_OUTPUT_STD  = Float32(2.4474483)

@propagate_inbounds function surface_roughness_ocean(ij, diagn, progn)
    surface = diagn.nlayers

    UVₛ = diagn.physics.surface_wind_speed[ij]
    Uₛ = diagn.grid.u_grid[ij, surface]
    Vₛ = diagn.grid.v_grid[ij, surface]
    ℵ = progn.ocean.sea_ice_concentration[ij]
    
    norm_ℵ   = (ℵ - O_INPUT_MEANS[1]) / O_INPUT_STDS[1]
    norm_Uₛ  = (Uₛ - O_INPUT_MEANS[2]) / O_INPUT_STDS[2]
    norm_Vₛ  = (Vₛ - O_INPUT_MEANS[3]) / O_INPUT_STDS[3]
    norm_UVₛ = (UVₛ - O_INPUT_MEANS[4]) / O_INPUT_STDS[4]

    predictors = [norm_ℵ  norm_Uₛ  norm_Vₛ  norm_UVₛ]
    input = Dict("x_input" => predictors)

    prediction = OCEAN_ROUGHNESS_MODEL_REF[](input)["linear_5"][1]
    log_surface_roughness = (prediction * O_OUTPUT_STD) + O_OUTPUT_MEAN
    surface_roughness = exp(log_surface_roughness)
    min_roughness = max(surface_roughness, 2.4361885e-5)
    return Float32(min_roughness)  # Min ERA5 value over ocean
end

@propagate_inbounds function surface_roughness_land(ij, diagn, progn)
    sd = progn.land.snow_depth[ij]
    g = diagn.dynamics.geopotential[ij]
    vₕ = diagn.physics.land.vegetation_high[ij]
    vₗ = diagn.physics.land.vegetation_low[ij]
    
    norm_vₕ = (vₕ - L_INPUT_MEANS[1]) / L_INPUT_STDS[1]
    norm_vₗ = (vₗ - L_INPUT_MEANS[2]) / L_INPUT_STDS[2]
    norm_g = (g - L_INPUT_MEANS[3]) / L_INPUT_STDS[3]
    norm_sd = (sd - L_INPUT_MEANS[4]) / L_INPUT_STDS[4]

    raw_predictors = [norm_vₕ norm_vₗ norm_g norm_sd]
    predictors = Float32.(real.(raw_predictors))
    input = Dict("float_input" => predictors)

    prediction = LAND_ROUGHNESS_MODEL_REF[](input)["variable"][1]
    log_surface_roughness = (prediction * L_OUTPUT_STD) + L_OUTPUT_MEAN
    surface_roughness = exp(log_surface_roughness)
    min_roughness = max(surface_roughness, 1.2999999e-3)
    return Float32(min_roughness)  # Min ERA5 value over land
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
