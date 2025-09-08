abstract type AbstractLandThermodynamics <: AbstractModelComponent end

export LandThermodynamics
@kwdef mutable struct LandThermodynamics{NF} <: AbstractLandThermodynamics
    "[OPTION] Heat conductivity λ_{soil} of the soil [W/(m K)]" 
    heat_conductivity::NF = 0.42

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

    "[OPTION] Soil wetness at wilting point [volume fraction]"
    wilting_point::NF = 0.17
	
    "[OPTION] Freezing temperature for snow on top soil layer [K]"
    snow_temp_freeze::NF = 273.15

    "[OPTION] Melting temperature for snow on top soil layer [K]"
    snow_temp_melt::NF = 273.15+2.5
	
    "[OPTION] Amount of cumulative snowfall required to trigger land snow scheme [m/m²]"
    snowfall_threshold::NF = 0.01
	
    "[OPTION] Fresh Snow albedo []"
    snow_albedo_fresh::NF = 0.8

    "[OPTION] Old Snow albedo []"
    snow_albedo_old::NF = 0.4
end

LandThermodynamics(SG::SpectralGrid; kwargs...) = LandThermodynamics{SG.NF}(; kwargs...)
initialize!(land::LandThermodynamics, model::PrimitiveEquation) = nothing
