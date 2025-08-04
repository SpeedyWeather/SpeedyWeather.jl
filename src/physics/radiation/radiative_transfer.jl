abstract type AbstractPlanckFunction end

@kwdef struct PlanckFunction{NF} <: AbstractPlanckFunction
    planck_constant::NF = 6.62607015e-34    # J/Hz = kg m^2 / s
    boltzmann_constant::NF = 1.380649e-23   # J/K
    speed_of_light::NF = 299792458.0        # m/s
    
    # precomputed
    hc::NF = planck_constant * speed_of_light
    hc²::NF = planck_constant * speed_of_light^2
end

PlanckFunction(SG::SpectralGrid; kwargs...) = SpectralGrid{SG.NF}(;kwargs...)

function (P::PlanckFunction)(T, nu)
    (; hc, hc²) = P
    k = P.boltzmann_constant

    return @. 2hc²*nu^3 / (exp(hc*nu / (k*T)) - 1)
end