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

# function barrier
parameterization!(ij, diagn, progn, sppt::StochasticallyPerturbedParameterizationTendencies, model) =
    sppt!(ij, diagn, sppt)

function sppt!(ij, diagn, sppt)
    
    r = diagn.grid.random_pattern[ij]
    (; taper) = sppt
    (; u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid) = diagn.tendencies

    @inbounds for k in eachlayer(u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid)
        R = 1 + r*taper[k]
        u_tend_grid[ij, k] *= R
        v_tend_grid[ij, k] *= R
        temp_tend_grid[ij, k] *= R
        humid_tend_grid[ij, k] *= R
    end
end