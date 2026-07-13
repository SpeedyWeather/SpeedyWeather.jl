# Flexible output of Terrarium state variables in coupled Terrarium-SpeedyWeather
# simulations. Constructor methods are added to the `SpeedyWeather.TerrariumOutput`
# function stub; they build `TerrariumOutputVariable`s with name, units, long name
# and dimensionality derived from Terrarium's variable descriptors.

import SpeedyWeather: TerrariumOutput, AbstractOutput, AbstractSimulation,
    path, path_or_nothing, output!, define_variable!, write_array!,
    hastime, is3D, NetCDFOutput
import SpeedyWeather.NCDatasets: NCDataset, defVar

# name of the vertical netCDF dimension shared by all 3D Terrarium output variables
const SOIL_DEPTH_DIM = "soil_depth"

"""Output variable for a variable of the Terrarium land state in
`simulation.variables.prognostic.land.terrarium`, used with the constructors
[`TerrariumOutput`](@ref). 3D (subsurface) variables are written on their own
vertical dimension `soil_depth` with the depths of the Terrarium soil layer
centres as coordinates; ocean points are filled with `missing_value`. Supported
by [`NetCDFOutput`](@ref), and by `ZarrOutput` once Zarr.jl is loaded.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct TerrariumOutputVariable{F} <: SpeedyWeather.AbstractOutputVariable
    "Name in the netCDF file, defaults to the Terrarium variable name"
    name::String
    "Units string, derived from the Terrarium variable's physical units"
    unit::String = "1"
    "Long name attribute, defaults to the Terrarium variable's description"
    long_name::String = name
    "Property path to the variable within the Terrarium state, e.g. `(:temperature,)`"
    path::Tuple{Vararg{Symbol}}
    "Which dimensions of (longitude, latitude, vertical, time) the variable has"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    "Value to fill non-land (ocean) grid points with"
    missing_value::Float64 = NaN
    "Lossless compression level of the netCDF variable"
    compression_level::Int = 3
    "Bitshuffle the data for compression"
    shuffle::Bool = true
    "Number of mantissa bits to keep for compression"
    keepbits::Int = 15
    "Element-wise transform applied before writing, e.g. unit conversion"
    transform::F = identity
    "Number of Terrarium soil layers (3D variables only)"
    nlayers::Int = 1
    "Depths of the soil layer centres [m], positive down, in field order (3D variables only)"
    depths::Vector{Float64} = Float64[]
    "Scratch field on the model grid to unmask the Terrarium land columns into"
    scratch_grid::Union{Nothing, RingGrids.AbstractField} = nothing
    "Scratch field on the output grid to interpolate into (3D variables only)"
    scratch_output::Union{Nothing, RingGrids.AbstractField} = nothing
end

# supported are 2D (XY) variables and 3D (XYZ) variables at cell centres;
# Face-located 3D variables (e.g. hydraulic_conductivity) have nlayers+1 points
supported_variable(descriptor) = supported_variable(Terrarium.vardims(descriptor))
supported_variable(::Terrarium.XY) = true
supported_variable(dims::Terrarium.XYZ) = dims.z isa Center

"""$(TYPEDSIGNATURES)
Create a [`TerrariumOutputVariable`](@ref) from a Terrarium variable `descriptor`
with name, units, long name and dimensionality derived from its metadata;
`kwargs` (e.g. `name`, `unit`, `keepbits`, `transform`) override the derived values."""
function terrarium_output_variable(
        model::Terrarium.AbstractModel,
        descriptor::Terrarium.AbstractProcessVariable;
        kwargs...
    )
    var_name = Terrarium.varname(descriptor)
    supported_variable(descriptor) ||
        throw(ArgumentError("Terrarium variable $var_name on $(Terrarium.vardims(descriptor)) " *
            "is not supported for output, only 2D (XY) and cell-centre 3D (XYZ) variables are."))

    units = Terrarium.varunits(descriptor)
    unit_string = string(units)
    unit = isempty(unit_string) ? "1" : unit_string     # NoUnits prints as ""
    long_name = isempty(descriptor.desc) ? replace(String(var_name), "_" => " ") : descriptor.desc

    dims = Terrarium.vardims(descriptor)
    is3d = dims isa Terrarium.XYZ
    field_grid = Terrarium.get_field_grid(model.grid)
    nlayers = is3d ? size(field_grid, 3) : 1
    # cell-centre depths [m], positive down; znodes is in field order (deepest layer first)
    depths = is3d ? Float64.(-znodes(field_grid, Center(), Center(), Center())) : Float64[]

    return TerrariumOutputVariable(;
        name = String(var_name), unit, long_name,
        path = (var_name,),
        dims_xyzt = (true, true, is3d, true),
        nlayers, depths,
        kwargs...
    )
end

# name-keyed NamedTuple of all top-level Terrarium variable descriptors of the
# selected groups; prognostic variables take precedence over auxiliary and input
# duplicates of the same name (as does `getproperty` on the Terrarium state)
function variable_descriptors(
        model::Terrarium.AbstractModel;
        prognostic::Bool = true,
        auxiliary::Bool = true,
        inputs::Bool = false,
    )
    tvars = Terrarium.Variables(Terrarium.variables(model))
    return merge(
        inputs ? tvars.inputs : (;),
        auxiliary ? tvars.auxiliary : (;),
        prognostic ? tvars.prognostic : (;),
    )
end

"""$(TYPEDSIGNATURES)
Create a [`TerrariumOutputVariable`](@ref) for the Terrarium variable `var_name`
of `model`, with units, long name and dimensionality derived from the variable's
metadata. Keyword arguments like `name` (netCDF name), `unit`, `keepbits` or
`transform` override the derived values."""
function SpeedyWeather.TerrariumOutput(
        model::Terrarium.AbstractModel,
        var_name::Symbol;
        kwargs...
    )
    descriptors = variable_descriptors(model, inputs = true)
    haskey(descriptors, var_name) ||
        throw(ArgumentError("Terrarium model has no variable $var_name, " *
            "available are: $(join(keys(descriptors), ", "))"))
    return terrarium_output_variable(model, descriptors[var_name]; kwargs...)
end

"""$(TYPEDSIGNATURES)
Create a tuple of [`TerrariumOutputVariable`](@ref)s for all (supported) variables
of the Terrarium `model`, to be added to a SpeedyWeather model's output via
`add!(model, TerrariumOutput(terrarium_model)...)`. The groups of Terrarium
variables to include are selected with `prognostic`, `auxiliary` and `inputs`
(inputs are excluded by default as they mirror the atmospheric forcing);
further keyword arguments are passed on to every single output variable."""
function SpeedyWeather.TerrariumOutput(
        model::Terrarium.AbstractModel;
        prognostic::Bool = true,
        auxiliary::Bool = true,
        inputs::Bool = false,
        kwargs...
    )
    descriptors = variable_descriptors(model; prognostic, auxiliary, inputs)
    supported = filter(supported_variable, collect(descriptors))
    return Tuple(terrarium_output_variable(model, descriptor; kwargs...) for descriptor in supported)
end

# convenience: construct directly from the SpeedyWeather land component
SpeedyWeather.TerrariumOutput(land::AbstractTerrariumLandModel, args...; kwargs...) =
    TerrariumOutput(land.model, args...; kwargs...)

# follow the property path into the Terrarium state (prognostic, auxiliary and
# input variables are all reachable via getproperty on StateVariables)
SpeedyWeather.path(variable::TerrariumOutputVariable, simulation) =
    foldl(getproperty, variable.path; init = simulation.variables.prognostic.land.terrarium)

"""$(TYPEDSIGNATURES)
Gather the Terrarium land columns of `variable` onto a full ring-grid scratch
field (ocean points remain at `missing_value`), interpolate onto the output
grid, transform and bitround. Returns `nothing` if `variable` is not part of
`simulation`, otherwise the output-grid field ready for a backend-specific
[`write_array!`](@ref) (shared between [`NetCDFOutput`](@ref) and `ZarrOutput`)."""
function terrarium_output_field!(
        output::AbstractOutput,
        variable::TerrariumOutputVariable,
        simulation::AbstractSimulation,
    )
    tfield = path_or_nothing(variable, simulation)      # Oceananigans field of the land columns
    isnothing(tfield) && return nothing                 # silently escape if not in simulation

    # lazily allocate the scratch fields on the first call
    output_NF = eltype(output.field2D)
    if isnothing(variable.scratch_grid)
        grid = SpeedyWeather.on_architecture(SpeedyWeather.CPU(), simulation.model.spectral_grid.grid)
        variable.scratch_grid = is3D(variable) ?
            RingGrids.Field(output_NF, grid, variable.nlayers) : RingGrids.Field(output_NF, grid)
        fill!(variable.scratch_grid, output_NF(variable.missing_value))
        if is3D(variable)
            variable.scratch_output = RingGrids.Field(output_NF, output.field2D.grid, variable.nlayers)
        end
    end

    # gather the land columns on CPU, drop the singleton y-dimension,
    # then scatter onto the (land points of the) full ring grid
    data = Array(interior(tfield))
    data = reshape(data, size(data, 1), size(data, 3))
    indices = SpeedyWeather.on_architecture(SpeedyWeather.CPU(), simulation.model.land.mask_indices)
    RingGrids.copy_unmasked!(variable.scratch_grid, data, indices)

    # interpolate onto the output grid, transform and bitround
    out = is3D(variable) ? variable.scratch_output : output.field2D
    RingGrids.interpolate!(out, variable.scratch_grid, output.interpolator)
    @. out = variable.transform(out)
    SpeedyWeather.round!(out, variable.keepbits)
    return out
end

"""$(TYPEDSIGNATURES)
Output `variable` into `output`: gather, interpolate, transform (via
[`terrarium_output_field!`](@ref)) and write to the netCDF file."""
function SpeedyWeather.output!(
        output::NetCDFOutput,
        variable::TerrariumOutputVariable,
        simulation::AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing
    out = terrarium_output_field!(output, variable, simulation)
    isnothing(out) && return nothing
    write_array!(output, variable, out)
    return nothing
end

# informative error for other output backends without a dedicated method
# (ZarrOutput is supported once Zarr.jl is loaded, see SpeedyWeatherTerrariumZarrExt)
SpeedyWeather.output!(output::AbstractOutput, variable::TerrariumOutputVariable, ::AbstractSimulation) =
    error("TerrariumOutput variables are not supported with $(nameof(typeof(output))). " *
        "Supported are NetCDFOutput, and ZarrOutput once Zarr.jl is loaded.")

"""$(TYPEDSIGNATURES)
Define `variable` in the netCDF `dataset`. 3D variables use their own vertical
dimension `soil_depth` (shared between all Terrarium output variables) with the
depths of the Terrarium soil layer centres [m, positive down] as coordinates."""
function SpeedyWeather.define_variable!(
        dataset::NCDataset,
        variable::TerrariumOutputVariable,
        output_NF::Type{<:AbstractFloat} = SpeedyWeather.DEFAULT_OUTPUT_NF,
    )
    if is3D(variable)   # lazily define the shared soil depth dimension + coordinates
        if !haskey(dataset.dim, SOIL_DEPTH_DIM)
            defVar(
                dataset, SOIL_DEPTH_DIM, variable.depths, (SOIL_DEPTH_DIM,),
                attrib = Dict(
                    "units" => "m", "long_name" => "depth of soil layer centre",
                    "positive" => "down"
                )
            )
        elseif dataset.dim[SOIL_DEPTH_DIM] != variable.nlayers
            error("Terrarium output variable $(variable.name) has $(variable.nlayers) soil layers, " *
                "but the $SOIL_DEPTH_DIM dimension already has $(dataset.dim[SOIL_DEPTH_DIM]).")
        end
    end

    attributes = Dict(
        "long_name" => variable.long_name, "units" => variable.unit,
        "_FillValue" => output_NF(variable.missing_value)
    )
    all_dims = ("lon", "lat", SOIL_DEPTH_DIM, "time")
    dims = collect(dim for (dim, this_dim) in zip(all_dims, variable.dims_xyzt) if this_dim)
    return defVar(
        dataset, variable.name, output_NF, dims, attrib = attributes;
        deflatelevel = variable.compression_level, shuffle = variable.shuffle
    )
end
