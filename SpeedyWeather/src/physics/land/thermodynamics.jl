abstract type AbstractLandThermodynamics <: AbstractModelComponent end

export LandThermodynamics
@kwdef mutable struct LandThermodynamics{NF} <: AbstractLandThermodynamics
    "[OPTION] Heat conductivity λ of the soil [W/(m K)]" 
    heat_conductivity::NF = 0.42

    "[OPTION] Field capacity (W_cap or γ) [volume fraction]"
    field_capacity::NF = 0.3

    "[OPTION] Heat capacity Cw of water [J/(m³ K)]"
    heat_capacity_water::NF = 4.2e6

    "[OPTION] Heat capacity Cs of dry soil [J/(m³ K)]"
    heat_capacity_dry_soil::NF = 1.13e6

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    wilting_point::NF = 0.17
end

LandThermodynamics(SG::SpectralGrid; kwargs...) = LandThermodynamics{SG.NF}(; kwargs...)
initialize!(land::LandThermodynamics, model::PrimitiveEquation) = nothing
