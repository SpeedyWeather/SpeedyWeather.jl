abstract type AbstractLandThermodynamics <: AbstractModelComponent end

export LandThermodynamics
@kwdef mutable struct LandThermodynamics{NF} <: AbstractLandThermodynamics
    "[OPTION] Heat conductivity λ of the soil [W/(m K)]" 
    heat_conductivity::NF = 0.42

    "[OPTION] Field capacity γ per meter soil [1]" #CPL changed this to be the same as below
    field_capacity::NF = 0.3

    "[OPTION] Heat capacity Cw of water [J/(m³ K)]"
    heat_capacity_water::NF = 4.2e6

    "[OPTION] Heat capacity Cs of dry soil [J/(m³ K)]"
    heat_capacity_dry_soil::NF = 1.13e6

    "[OPTION] Soil wetness at field capacity [volume fraction]"
    W_cap::NF = 0.3

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    W_wilt::NF = 0.17
end

LandThermodynamics(SG::SpectralGrid; kwargs...) = LandThermodynamics{SG.NF}(; kwargs...)
initialize!(land::LandThermodynamics, model::PrimitiveEquation) = nothing
