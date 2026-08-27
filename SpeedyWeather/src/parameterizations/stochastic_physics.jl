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
    taper = on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers))
    return StochasticallyPerturbedParameterizationTendencies(tapering, taper)
end

Adapt.@adapt_structure StochasticallyPerturbedParameterizationTendencies

# No additional variables required, return empty tuple
variables(::AbstractStochasticPhysics) = ()

# precompute the taper
function initialize!(sppt::StochasticallyPerturbedParameterizationTendencies, model::PrimitiveEquation)
    # SPPT perturbs with vars.grid.random_pattern but doesn't define that variable itself, it has
    # to come from another component, e.g. a random process. Check here, otherwise sppt! errors
    # only later inside the column parameterizations (on GPU as an unsupported jl_f_getfield call)
    has_random_pattern = any(all_variables(model)) do var
        var isa GridVariable && var.name == :random_pattern && var.namespace == Symbol()
    end
    @assert has_random_pattern "StochasticallyPerturbedParameterizationTendencies requires a "*
        "`random_pattern` grid variable, define a random process for the model, e.g. "*
        "`random_process = SpectralAR1Process(spectral_grid)`."

    coord = model.geometry.vertical_coordinates
    nlayers = get_nlayers(coord)
    (; taper) = sppt
    arch = architecture(taper)

    launch!(
        arch, LinearWorkOrder, (nlayers,), initialize_sppt_taper_kernel!,
        taper, sppt.tapering, coord
    )
    return nothing
end

@kernel inbounds = true function initialize_sppt_taper_kernel!(taper, tapering, coordinate)
    k = @index(Global, Linear)
    taper[k] = tapering(sigma(k, coordinate))
end

# function barrier
@propagate_inbounds parameterization!(ij, vars, sppt::StochasticallyPerturbedParameterizationTendencies, model) =
    sppt!(ij, vars, sppt, model.time_stepping)

"""$(TYPEDSIGNATURES)
Apply stochastically perturbed parameterization tendencies (SPPT) to
u, v, temperature and humidity in column ij."""
@propagate_inbounds function sppt!(ij, vars, sppt, time_stepping)

    r = vars.grid.random_pattern[ij]
    (; taper) = sppt
    u_tend = get_tendency_step(vars.tendencies.grid.u, time_stepping, sppt)
    v_tend = get_tendency_step(vars.tendencies.grid.v, time_stepping, sppt)
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, sppt)

    # dry models don't have humidity just perturb a dummy array to avoid branching in the loop below
    humid_tend = haskey(vars.tendencies.grid, :humidity) ?
        get_tendency_step(vars.tendencies.grid.humidity, time_stepping, sppt) : vars.scratch.a_grid


    @inbounds for k in eachlayer(u_tend, v_tend, temp_tend, humid_tend)
        R = 1 + r * taper[k]        # r in [-1, 1], R in [0, 2] (don't change sign of tendency)
        u_tend[ij, k] *= R          # perturb all prognostic variables in the same way
        v_tend[ij, k] *= R
        temp_tend[ij, k] *= R
        humid_tend[ij, k] *= R
    end
    return nothing
end
