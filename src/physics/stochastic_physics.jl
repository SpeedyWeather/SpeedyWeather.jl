abstract type AbstractStochasticPhysics <: AbstractParameterization end

# fucntion barriers
function perturb_parameterization_inputs!(column::AbstractColumnVariables, model::PrimitiveEquation)
    perturb_parameterization_inputs!(column, model.stochastic_physics, model)
end

function perturb_parameterization_tendencies!(column::AbstractColumnVariables, model::PrimitiveEquation)
    perturb_parameterization_tendencies!(column, model.stochastic_physics, model)
end


export NoStochasticPhysics
struct NoStochasticPhysics <: AbstractStochasticPhysics end
NoStochasticPhysics(::SpectralGrid) = NoStochasticPhysics()
initialize!(::NoStochasticPhysics, ::PrimitiveEquation) = nothing
perturb_parameterization_inputs!(::ColumnVariables, ::NoStochasticPhysics, ::PrimitiveEquation) = nothing
perturb_parameterization_tendencies!(::ColumnVariables, ::NoStochasticPhysics, ::PrimitiveEquation) = nothing

export StochasticallyPerturbedParameterizationTendencies
@kwdef struct StochasticallyPerturbedParameterizationTendencies{NF, VectorType} <: AbstractStochasticPhysics
    
    "Number of vertical layers"
    nlayers::Int

    "[OPTION] Vertical tapering function, reduce strength towards surface (σ=1)"
    tapering::Function = σ -> 1 # σ < 0.8 ? 1 : 1 - (σ - 0.8)/0.2

    "[DERIVED] Precalculate vertical tapering during initialization"
    taper::VectorType = zeros(NF, nlayers)
end

# generator function
function StochasticallyPerturbedParameterizationTendencies(SG::SpectralGrid; kwargs...)
    (; nlayers, NF, VectorType) = SG
    return StochasticallyPerturbedParameterizationTendencies{NF, VectorType}(; nlayers, kwargs...)
end

function initialize!(sppt::StochasticallyPerturbedParameterizationTendencies, model::PrimitiveEquation)
    sppt.taper .= sppt.tapering.(model.geometry.σ_levels_full)
    return nothing
end

# only perturb tendencies (=outputs) not inputs
function perturb_parameterization_inputs!(
    ::AbstractColumnVariables,
    ::StochasticallyPerturbedParameterizationTendencies,
    ::PrimitiveEquation)
    return nothing
end

function perturb_parameterization_tendencies!(
    column::AbstractColumnVariables,
    sppt::StochasticallyPerturbedParameterizationTendencies,
    model::PrimitiveEquation)
    
    r = column.random_value
    (; taper) = sppt
    (; u_tend, v_tend, temp_tend, humid_tend) = column

    for k in eachlayer(column)
        R = 1 + r*taper[k]
        u_tend[k] *= R
        v_tend[k] *= R
        temp_tend[k] *= R
        humid_tend[k] *= R
    end
end