names = (
    (:OutgoingShortwaveRadiationOutput, "osr",  "Outgoing shortwave radiation",     :outgoing_shortwave_radiation),
    (:OutgoingLongwaveRadiationOutput,  "olr",  "Outgoing longwave radiation",      :outgoing_longwave_radiation),
    (:SurfaceShortwaveUpOutput,         "sru",  "Surface shortwave radiation up",   :surface_shortwave_up),
    (:SurfaceShortwaveDownOutput,       "srd",  "Surface shortwave radiation down", :surface_shortwave_down),
    (:SurfaceLongwaveUpOutput,          "lru",  "Surface longwave radiation up",    :surface_longwave_up),
    (:SurfaceLongwaveDownOutput,        "lrd",  "Surface longwave radiation down",  :surface_longwave_down),
)

for name in names
    typename, shortname, longname, varname = name
    @eval begin
        @kwdef mutable struct $typename <: AbstractOutputVariable
            name::String = $shortname
            unit::String = "W/m^2"
            long_name::String = $longname
            dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
            missing_value::Float64 = NaN
            compression_level::Int = 3
            shuffle::Bool = true
            keepbits::Int = 7
        end

        path(::$typename, simulation) =
            simulation.diagnostic_variables.physics.$varname
    end
end

RadiationOutput() = (
    OutgoingLongwaveRadiationOutput(),
    OutgoingShortwaveRadiationOutput(),
    SurfaceShortwaveUpOutput(),
    SurfaceShortwaveDownOutput(),
    SurfaceLongwaveUpOutput(),
    SurfaceLongwaveDownOutput(),
)