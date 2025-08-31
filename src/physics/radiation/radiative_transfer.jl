abstract type AbstractPlanckFunction <: AbstractModelComponent end

"""$(TYPEDSIGNATURES)
Radiance as function of wavenumber (1/m) and temperature (K) following
Planck's law."""
@kwdef struct PlanckFunction{NF} <: AbstractPlanckFunction
    "Planck constant [J/Hz]"
    planck_constant::NF = 6.62607015e-34    # J/Hz = kg m^2 / s

    "Boltzmann constant [J/K]"
    boltzmann_constant::NF = 1.380649e-23   # J/K

    "Speed of light [m/s]"
    speed_of_light::NF = 299792458.0        # m/s
    
    "[DERIVED] Planck constant times speed of light squared"
    hc²::NF = planck_constant * speed_of_light^2
    
    "[DERIVED] Planck constant times speed of light"
    hc::NF = planck_constant * speed_of_light

    "[DERIVED] hc divided by Boltzmann constant"
    hc_kB::NF = hc / boltzmann_constant
end

PlanckFunction(SG::SpectralGrid; kwargs...) = PlanckFunction{SG.NF}(; kwargs...)

struct PlanckWavelength end
struct PlanckWavenumber end

(P::PlanckFunction)(nu, T) = P(nu, T, PlanckWavenumber)

"""$(TYPEDSIGNATURES)
Planck function from wavenumber `nu` (1/m) and temperature `T` (K),
applied as functor from `P::PlanckFunction`."""
function (P::PlanckFunction)(nu, T, ::Type{PlanckWavenumber})
    (; hc_kB, hc²) = P
    return @. 2hc²*nu^3 / (exp(hc_kB*(nu/T)) - 1)
end

function (P::PlanckFunction)(λ, T, ::Type{PlanckWavelength})
    (; hc_kB, hc²) = P
    return @. 2hc² / (λ^5*(exp(hc_kB/(λ*T)) - 1))
end