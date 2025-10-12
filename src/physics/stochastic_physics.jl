abstract type AbstractStochasticPhysics <: AbstractParameterization end

# fucntion barriers
function perturb_inputs!(ij, diagn, progn, model)
    perturb_inputs!(ij, diagn, progn, model.stochastic_physics, model)
end

function perturb_tendencies!(ij, diagn, progn, model)
    perturb_tendencies!(ij, diagn, progn, model.stochastic_physics, model)
end

# no perturbations
perturb_inputs!(ij, diagn, progn, ::Nothing, model) = nothing
perturb_tendencies!(ij, diagn, progn, ::Nothing, model) = nothing

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
function StochasticallyPerturbedParameterizationTendencies(SG::SpectralGrid; 
    tapering = σ -> 1, # σ < 0.8 ? 1 : 1 - (σ - 0.8)/0.2
)
    taper = on_architecture(SG.architecture, zeros(SG.nlayers))
    return StochasticallyPerturbedParameterizationTendencies(tapering, taper)
end

function initialize!(sppt::StochasticallyPerturbedParameterizationTendencies, model::PrimitiveEquation)
    sppt.taper .= sppt.tapering.(model.geometry.σ_levels_full)
    return nothing
end

# only perturb tendencies (=outputs) not inputs
perturb_inputs!(ij, diagn, progn, sppt::StochasticallyPerturbedParameterizationTendencies, model) = nothing

function perturb_tendencies!(ij, diagn, progn, sppt::StochasticallyPerturbedParameterizationTendencies, model)
    
    r = diagn.grid.random_pattern[ij]
    (; taper) = sppt
    (; u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid) = diagn.tendencies

    for k in eachlayer(u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid)
        R = 1 + r*taper[k]
        u_tend_grid[ij, k] *= R
        v_tend_grid[ij, k] *= R
        temp_tend_grid[ij, k] *= R
        humid_tend_grid[ij, k] *= R
    end
end