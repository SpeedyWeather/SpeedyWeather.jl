names = (
    (:OutgoingShortwaveRadiation, "osr", "Outgoing shortwave radiation", :outgoing_shortwave),
    (:OutgoingLongwaveRadiation, "olr", "Outgoing longwave radiation", :outgoing_longwave),
    (:SurfaceShortwaveUp, "sru", "Surface shortwave radiation up", :surface_shortwave_up),
    (:SurfaceShortwaveDown, "srd", "Surface shortwave radiation down", :surface_shortwave_down),
    (:SurfaceLongwaveUp, "lru", "Surface longwave radiation up", :surface_longwave_up),
    (:SurfaceLongwaveDown, "lrd", "Surface longwave radiation down", :surface_longwave_down),
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
            compression_level::Int = DEFAULT_COMPRESSION_LEVEL
            shuffle::Bool = DEFAULT_SHUFFLE
            keepbits::Int = 7
        end

        path(::$typename, simulation) =
            simulation.variables.parameterizations.$varname
    end
end

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct Albedo <: AbstractOutputVariable
    name::String = "albedo"
    unit::String = "1"
    long_name::String = "albedo"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)  # albedo is diagnostic
    missing_value::Float64 = NaN
    compression_level::Int = DEFAULT_COMPRESSION_LEVEL
    shuffle::Bool = DEFAULT_SHUFFLE
    keepbits::Int = DEFAULT_KEEPBITS
end

path(::Albedo, simulation) = simulation.variables.parameterizations.albedo

Radiation() = (
    OutgoingLongwaveRadiation(),
    OutgoingShortwaveRadiation(),
    SurfaceShortwaveUp(),
    SurfaceShortwaveDown(),
    SurfaceLongwaveUp(),
    SurfaceLongwaveDown(),
    Albedo(),
)
