export ModelLayers, PressureLayers

"""Supertype of the vertical layers that an output writer writes its 3D variables on.
A property of the writer (like its output grid), not of the individual output variables,
set with the `layers` keyword argument, e.g. `NetCDFOutput(spectral_grid, layers =
PressureLayers(spectral_grid))`. Subtypes are [`ModelLayers`](@ref) and
[`PressureLayers`](@ref)."""
abstract type AbstractOutputLayers <: AbstractModelComponent end

"""Write 3D variables on the model's vertical layers, sigma or hybrid sigma-pressure
depending on `model.geometry.vertical_coordinates`. Default for all output writers."""
struct ModelLayers <: AbstractOutputLayers end

"""Pressure [Pa] of the layers that `PressureLayers` interpolates onto by default."""
const DEFAULT_PRESSURE_LAYERS = [50, 100, 200, 300, 500, 700, 850, 925, 1000] .* 100.0

"""Write 3D atmospheric variables interpolated onto `pressure` layers [Pa] instead of the
model's vertical layers, see [Vertical interpolation onto pressure layers](@ref
vertical_interpolation). Construct with a `SpectralGrid` and pass to an output writer as
`NetCDFOutput(spectral_grid, layers = PressureLayers(spectral_grid))`. Land (soil)
variables and 2D variables are unaffected. Fields are $(TYPEDFIELDS)"""
struct PressureLayers{V, I, F} <: AbstractOutputLayers
    "[OPTION] pressure layers [Pa] to interpolate onto, monotonic, e.g. ascending from the
    top of the atmosphere to the surface. In the model's number format and on its
    architecture so that the interpolation kernel can run on the GPU"
    pressure::V

    "[OPTION] interpolate linearly in pressure or in its logarithm?"
    interpolation::I

    "[DERIVED] scratch field (npoints, npressure) on the model grid to interpolate into
    before the horizontal interpolation onto the output grid"
    scratch::F
end

"""$(TYPEDSIGNATURES)
Vertical output layers at the given `pressure` [Pa], $(DEFAULT_PRESSURE_LAYERS ./ 100) hPa
by default. Allocates the pressure vector and one scratch field (shared by all output
variables) on `SG`'s architecture and in its number format, so that the vertical
interpolation runs where the model runs."""
function PressureLayers(
        SG::SpectralGrid;
        pressure::AbstractVector = DEFAULT_PRESSURE_LAYERS,
        interpolation::AbstractVerticalInterpolation = LinearInLogPressure(),
    )
    # the interpolation itself doesn't care about the order (every layer is searched for
    # independently) but a non-monotonic vertical coordinate in the output file would be
    # a mistake, so require monotonic, ascending (top to surface) or descending
    monotonic = issorted(pressure) || issorted(pressure, rev = true)
    @assert monotonic "Output pressure layers have to be monotonic, got $pressure"
    @assert all(>(0), pressure) "Output pressure layers have to be positive, got $pressure"

    (; NF, grid) = SG
    return PressureLayers(
        on_architecture(SG.architecture, NF.(pressure)),
        interpolation,
        Field(NF, grid, length(pressure)),      # on the model grid, not the output grid
    )
end

# convenience: PressureLayers(spectral_grid, [850, 500, 200] .* 100)
PressureLayers(SG::SpectralGrid, pressure::AbstractVector; kwargs...) =
    PressureLayers(SG; pressure, kwargs...)

function Base.show(io::IO, layers::PressureLayers)
    pressure = on_architecture(CPU(), layers.pressure) ./ 100
    println(io, styled"{warning:PressureLayers}")
    println(io, styled"├ {info:pressure} = $pressure hPa")
    return print(io, styled"└ {info:interpolation}::$(typeof(layers.interpolation))")
end

Base.show(io::IO, ::ModelLayers) = print(io, styled"{warning:ModelLayers}")

"""$(TYPEDSIGNATURES)
Number of vertical layers that 3D atmospheric variables are written on, used by the output
writers to allocate their 3D scratch field on the output grid."""
get_nlayers(::ModelLayers, SG::SpectralGrid) = SG.nlayers
get_nlayers(layers::PressureLayers, ::SpectralGrid) = length(layers.pressure)

"""$(TYPEDSIGNATURES)
Name of the vertical dimension in the output file or store for 3D atmospheric variables."""
vertical_dimension_name(::ModelLayers) = "layer"
vertical_dimension_name(::PressureLayers) = "pressure"

"""$(TYPEDSIGNATURES)
Vertical dimension of `variable` in `output`. Variables on the model's atmospheric layers
follow the output writer's `layers`, e.g. "layer" or "pressure". Variables with their own
vertical dimension (soil variables, custom output variables that extend the 1-argument
[`vertical_dimension_name`](@ref)) keep theirs."""
function vertical_dimension_name(output::AbstractOutput, variable::AbstractOutputVariable)
    dim = vertical_dimension_name(variable)
    return dim == "layer" ? vertical_dimension_name(output_layers(output)) : dim
end

"""$(TYPEDSIGNATURES)
Vertical output layers of `output`, [`ModelLayers`](@ref) for output writers that don't
define a `layers` field."""
output_layers(output::AbstractOutput) = hasproperty(output, :layers) ? output.layers : ModelLayers()

"""$(TYPEDSIGNATURES)
Extrapolation of `variable` beyond the model's vertical layers: its `extrapolation` field
if it has one (e.g. [`TemperatureOutput`](@ref) descends dry-adiabatically below the
lowest model layer), otherwise constant. Same `hasproperty` pattern as `transform` and
`unscale` in the generic [`output!`](@ref)."""
output_extrapolation(variable::AbstractOutputVariable) =
    hasproperty(variable, :extrapolation) ? variable.extrapolation : ConstantExtrapolation()

"""$(TYPEDSIGNATURES)
Interpolate `field` onto the output writer's vertical `layers`, returning the field to
write out. No-op on [`ModelLayers`](@ref) and for variables that are not 3D atmospheric
variables (2D variables, soil variables). For [`PressureLayers`](@ref) this interpolates
into the layers' scratch field, on the model grid and the model's architecture, i.e.
before the horizontal interpolation onto the output grid."""
interpolate_layers!(::ModelLayers, field, ::AbstractOutputVariable, ::AbstractSimulation) = field

function interpolate_layers!(
        layers::PressureLayers,
        field::AbstractField,
        variable::AbstractOutputVariable,
        simulation::AbstractSimulation,
    )
    # only 3D atmospheric variables, soil variables keep their own vertical dimension
    (is3D(variable) && !is_land(variable)) || return field

    (; model) = simulation
    surface_pressure = simulation.variables.parameterizations.surface_pressure

    interpolate_pressure_layers!(
        layers.scratch, field, surface_pressure, layers.pressure,
        model.geometry.vertical_coordinates, layers.interpolation,
        output_extrapolation(variable),
    )
    return layers.scratch
end

"""$(TYPEDSIGNATURES)
Set κ of any [`DryAdiabaticExtrapolation`](@ref) in `output.variables` from the model's
atmosphere, so that the adiabatic descent below the lowest model layer uses the same
κ = R_dry/cₚ as the model itself. Called at initialize!, no-op for a model without an
atmosphere (Barotropic) or an output writer that has no output variables (JLD2Output)."""
function sync_extrapolations!(output::AbstractOutput, model::AbstractModel)
    hasproperty(model, :atmosphere) && hasproperty(output, :variables) || return nothing
    for variable in values(output.variables)
        hasproperty(variable, :extrapolation) || continue
        variable.extrapolation = sync_extrapolation(variable.extrapolation, model.atmosphere)
    end
    return nothing
end

sync_extrapolation(E::AbstractVerticalExtrapolation, ::AbstractAtmosphere) = E

# keep the number format of the extrapolation so that the type of `variable.extrapolation`
# doesn't change, κ is converted to the field's number format in `extrapolate_below` anyway
sync_extrapolation(::DryAdiabaticExtrapolation{NF}, atmosphere::AbstractAtmosphere) where {NF} =
    DryAdiabaticExtrapolation(NF(atmosphere.κ))
sync_extrapolation(E::SubsurfaceMask, atmosphere::AbstractAtmosphere) =
    SubsurfaceMask(sync_extrapolation(E.above_surface, atmosphere), E.missing_value)

"""$(TYPEDSIGNATURES)
Define the vertical coordinate of 3D atmospheric variables in the output file or store
`dest`: sigma for [`ModelLayers`](@ref), pressure in hPa for [`PressureLayers`](@ref).
Backend-agnostic via [`define_coordinate!`](@ref)."""
function define_vertical_coordinate!(dest, layers::ModelLayers, model::AbstractModel)
    σ = on_architecture(CPU(), model.geometry.σ_levels_full)
    return define_coordinate!(
        dest, vertical_dimension_name(layers), collect(Float64.(σ)),
        attribs = Dict("units" => "1", "long_name" => "sigma layer"),
    )
end

function define_vertical_coordinate!(dest, layers::PressureLayers, ::AbstractModel)
    pressure = on_architecture(CPU(), layers.pressure)       # coordinates are written from CPU
    return define_coordinate!(
        dest, vertical_dimension_name(layers), collect(Float64.(pressure) ./ 100),
        attribs = Dict(
            "units" => "hPa", "long_name" => "pressure", "positive" => "down",
            "standard_name" => "air_pressure",
        ),
    )
end
