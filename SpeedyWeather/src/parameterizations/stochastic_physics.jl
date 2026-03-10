abstract type AbstractStochasticPhysics <: AbstractParameterization end

export StochasticallyPerturbedParameterizationTendencies

"""Defines the stochastically perturbed parameterization tendencies (SPPT)
including an optional tapering as a function of the vertical sigma level
$(TYPEDFIELDS)"""
struct StochasticallyPerturbedParameterizationTendencies{F, VectorType} <: AbstractStochasticPhysics
    "[OPTION] Vertical tapering function, reduce strength towards surface (σ=1)"
    tapering::F

    "[DERIVED] Precalculate vertical tapering during initialization"
    taper::VectorType
end

# generator function
function StochasticallyPerturbedParameterizationTendencies(
        SG::SpectralGrid;
        tapering = σ -> 1, # σ < 0.8 ? 1 : 1 - (σ - 0.8)/0.2
    )
    taper = on_architecture(SG.architecture, zeros(SG.nlayers))
    return StochasticallyPerturbedParameterizationTendencies(tapering, taper)
end

Adapt.@adapt_structure StochasticallyPerturbedParameterizationTendencies

# No additional variables required, return empty tuple
variables(::AbstractStochasticPhysics) = ()

# precompute the taper
function initialize!(sppt::StochasticallyPerturbedParameterizationTendencies, model::PrimitiveEquation)
    sppt.taper .= sppt.tapering.(model.geometry.σ_levels_full)
    return nothing
end

# function barrier
@propagate_inbounds parameterization!(ij, vars, sppt::StochasticallyPerturbedParameterizationTendencies, model) =
    sppt!(ij, vars, sppt)

"""$(TYPEDSIGNATURES)
Apply stochastically perturbed parameterization tendencies (SPPT) to
u, v, temperature and humidity in column ij."""
@propagate_inbounds function sppt!(ij, vars, sppt)

    r = vars.grid.random_pattern[ij]
    (; taper) = sppt
    u_tend = vars.tendencies.grid.u
    v_tend = vars.tendencies.grid.v
    temp_tend = vars.tendencies.grid.temp

    # dry models don't have humidity just perturb a dummy array to avoid branching in the loop below
    humid_tend = haskey(vars.tendencies.grid, :humid) ? vars.tendencies.grid.humid : vars.scratch.a_grid


    @inbounds for k in eachlayer(u_tend, v_tend, temp_tend, humid_tend)
        R = 1 + r * taper[k]        # r in [-1, 1], R in [0, 2] (don't change sign of tendency)
        u_tend[ij, k] *= R          # perturb all prognostic variables in the same way
        v_tend[ij, k] *= R
        temp_tend[ij, k] *= R
        humid_tend[ij, k] *= R
    end
    return nothing
end
