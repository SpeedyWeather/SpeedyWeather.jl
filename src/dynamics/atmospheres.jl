"""
    EA = EarthAtmosphere(kwargs...)

Create a struct `EarthAtmosphere<:AbstractPlanet`, with the following physical/chemical
characteristics. Note that `radius` is not part of it as this should be chosen
in `SpectralGrid`. Default keyword arguments are

$(TYPEDFIELDS)"""
@kwdef struct EarthAtmosphere <: AbstractAtmosphere
    # ATMOSPHERE

    "molar mass of dry air [g/mol]"
    mol_mass_dry_air::Float64 = 28.9649

    "molar mass of water vapour [g/mol]"
    mol_mass_vapour::Float64 = 18.0153

    "specific heat at constant pressure [J/K/kg]"
    cₚ::Float64 = 1004

    "universal gas constant [J/K/mol]"
    R_gas::Float64 = 8.3145

    "specific gas constant for dry air [J/kg/K]"
    R_dry::Float64 = 1000*R_gas/mol_mass_dry_air

    "specific gas constant for water vapour [J/kg/K]"
    R_vapour::Float64 = 1000*R_gas/mol_mass_vapour

    "latent heat of condensation [J/g] for consistency with specific humidity [g/Kg]"
    alhc::Float64 = 2501

    "latent heat of sublimation [?]"
    alhs::Float64 = 2801

    "stefan-Boltzmann constant [W/m²/K⁴]"
    sbc::Float64 = 5.67e-8


    # STANDARD ATMOSPHERE (reference values)

    "moist adiabatic temperature lapse rate ``-dT/dz`` [K/km]"
    lapse_rate::Float64 = 5

    "absolute temperature at surface ``z=0`` [K]"
    temp_ref::Float64 = 288

    "absolute temperature in stratosphere [K]"
    temp_top::Float64 = 216

    "for stratospheric lapse rate [K] after Jablonowski"
    ΔT_stratosphere::Float64 = 4.8e5

    "start of the stratosphere in sigma coordinates"
    σ_tropopause::Float64 = 0.2

    "scale height for pressure [km]"
    scale_height::Float64 = 7.5

    "surface pressure [hPa]"
    pres_ref::Float64 = 1000

    "scale height for specific humidity [km]"
    scale_height_humid::Float64 = 2.5

    "relative humidity of near-surface air [1]"
    relhumid_ref::Float64 = 0.7

    "saturation water vapour pressure [Pa]"
    water_pres_ref::Float64 = 17

    "layer thickness for the shallow water model [km]"
    layer_thickness::Float64 = 8.5
end

function Base.show(io::IO,atm::AbstractAtmosphere)
    println(io,"$(typeof(atm))(")
    fields = propertynames(atm)
    nfields = length(fields)
    for i in 1:nfields
        key = fields[i]
        val = getfield(atm,key)
        s = "  $key::$(typeof(val)) = $val"
        if i < nfields println(io,s) else print(io,s*")") end
    end
end