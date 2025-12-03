abstract type AbstractLandThermodynamics <: AbstractModelComponent end

export LandThermodynamics
@kwdef mutable struct LandThermodynamics{NF} <: AbstractLandThermodynamics
    "[OPTION] Heat conductivity λ_{soil} of the soil [W/(m K)]" 
    heat_conductivity_dry_soil::NF = 0.42

    "[OPTION] Heat conductivity λ_{snow} of snow [W/(m K)]" 
    heat_conductivity_snow::NF = 0.08

    "[OPTION] Field capacity (W_cap or γ) [volume fraction]"
    field_capacity::NF = 0.3

    "[OPTION] Heat capacity Cw of water [J/(m³ K)]"
    heat_capacity_water::NF = 4.2e6

    "[OPTION] Heat capacity Ci of snow [J/(m³ K)]"
    heat_capacity_snow::NF = 2.09e3

    "[OPTION] Heat capacity Cs of dry soil [J/(m³ K)]"
    heat_capacity_dry_soil::NF = 1.13e6

    "[OPTION] Soil density [kg/m³]"
    soil_density::NF = 1700

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    wilting_point::NF = 0.17
end

LandThermodynamics(SG::SpectralGrid, geometry::LandGeometry; kwargs...) = LandThermodynamics{SG.NF}(; kwargs...)
initialize!(land::LandThermodynamics, model::PrimitiveEquation) = nothing
