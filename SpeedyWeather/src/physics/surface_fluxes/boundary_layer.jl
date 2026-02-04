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

export LearnedSurfaceRoughness
@kwdef struct LearnedSurfaceRoughness{NF, LM, LNN, LP, LS} <: AbstractSurfaceRoughness
    # Ocean normalisation parameters
    ocean_input_means::Vector{NF} = Float32[0.19490805, 0.11980359, 7.7569385]
    ocean_input_stds::Vector{NF} = Float32[6.6221075, 5.4018555, 3.5962458]
    ocean_output_mean::NF = -9.1571865f0
    ocean_output_std::NF = 1.0090618f0

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
    land_output_mean::NF = -5.031811f0
    land_output_std::NF = 2.4447718f0

    land_input_buffer::Vector{Float32} = zeros(Float32, 13)

    land_model::LM
    land_nn::LNN
    land_params::LP
    land_states::LS
end

function Base.show(io::IO, scheme::LearnedSurfaceRoughness)
    print(io, "LearnedSurfaceRoughness{$(eltype(scheme.ocean_output_mean))}")
    println(io)
    println(io, "├ Ocean: Empirically derived analytical model")
    n_layers = length(keys(scheme.land_params))
    println(io, "└ Land:  Neural network ($n_layers layers)")
end

Adapt.@adapt_structure ConstantSurfaceRoughness
ConstantSurfaceRoughness(SG::SpectralGrid, kwargs...) = ConstantSurfaceRoughness{SG.NF}(; kwargs...)
initialize!(::ConstantSurfaceRoughness, ::PrimitiveEquation) = nothing

@propagate_inbounds function surface_roughness_land(ij, scheme::ConstantSurfaceRoughness, diagn, progn)
    return scheme.roughness_length_land
end

@propagate_inbounds function surface_roughness_ocean(ij, scheme::ConstantSurfaceRoughness, diagn, progn)
    return scheme.roughness_length_ocean
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

function Base.show(io::IO, B::BulkRichardsonDrag)
    println(io, "BulkRichardsonDrag{$(typeof(B.von_Karman))}")
    println(io, "├ von_Karman: $(B.von_Karman)")
    println(io, "├ critical_Richardson: $(B.critical_Richardson)")
    println(io, "├ drag_min: $(B.drag_min)")
    print(io,   "└┐surface_roughness: ")
    
    buf = IOBuffer()
    show(buf, B.surface_roughness)
    s = String(take!(buf))
    lines = split(s, '\n')
    
    if !isempty(lines)
        println(io, lines[1])
    end

    prefix = " " 
    for line in lines[2:end]
        if !isempty(line)
            println(io, prefix, line)
        end
    end
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
    z₀_land = zero(land_fraction)
    z₀_ocean = zero(land_fraction)
    z₀_land = land_fraction > 0 ? surface_roughness_land(ij, drag.surface_roughness, diagn, progn) : zero(land_fraction)
    z₀_ocean = land_fraction < 1 ? surface_roughness_ocean(ij, drag.surface_roughness, diagn, progn) : zero(land_fraction)

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
