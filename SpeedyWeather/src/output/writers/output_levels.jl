export ModelLevels, PressureLevels

"""Supertype of the vertical levels that an output writer writes its 3D variables on.
A property of the writer (like its output grid), not of the individual output variables,
set with the `levels` keyword argument, e.g. `NetCDFOutput(spectral_grid, levels =
PressureLevels())`. Subtypes are [`ModelLevels`](@ref) and [`PressureLevels`](@ref)."""
abstract type AbstractOutputLevels <: AbstractModelComponent end

"""Write 3D variables on the model's vertical levels, sigma or hybrid sigma-pressure
depending on `model.geometry.vertical_coordinates`. Default for all output writers."""
struct ModelLevels <: AbstractOutputLevels end

"""Write 3D atmospheric variables interpolated onto `pressure` levels [Pa] instead of the
model's vertical levels, see [Vertical interpolation onto pressure levels](@ref
vertical_interpolation). Pass to an output writer as
`NetCDFOutput(spectral_grid, levels = PressureLevels())`. Land (soil) variables and 2D
variables are unaffected. Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct PressureLevels{V, I} <: AbstractOutputLevels
    "[OPTION] pressure levels [Pa] to interpolate onto, monotonic, e.g. ascending from the top of the atmosphere to the surface"
    pressure::V = [50, 100, 200, 300, 500, 700, 850, 925, 1000] .* 100.0

    "[OPTION] interpolate linearly in pressure or in its logarithm?"
    interpolation::I = LinearInLogPressure()

    "[DERIVED] pressure levels in the model's number format and on its architecture,
    allocated at initialize!"
    pressure_on_arch::Union{AbstractVector, Nothing} = nothing

    "[DERIVED] scratch field (npoints, npressure) on the model grid to interpolate into
    before the horizontal interpolation onto the output grid, allocated at initialize!"
    scratch::Union{AbstractField, Nothing} = nothing
end

# convenience: PressureLevels([850, 500, 200] .* 100)
PressureLevels(pressure::AbstractVector; kwargs...) = PressureLevels(; pressure, kwargs...)

function Base.show(io::IO, levels::PressureLevels)
    println(io, styled"{warning:PressureLevels}")
    println(io, styled"├ {info:pressure} = $(levels.pressure ./ 100) hPa")
    return print(io, styled"└ {info:interpolation}::$(typeof(levels.interpolation))")
end

Base.show(io::IO, ::ModelLevels) = print(io, styled"{warning:ModelLevels}")

"""$(TYPEDSIGNATURES)
Number of vertical levels that 3D atmospheric variables are written on, used by the output
writers to allocate their 3D scratch field on the output grid."""
get_nlayers(::ModelLevels, SG::SpectralGrid) = SG.nlayers
get_nlayers(levels::PressureLevels, ::SpectralGrid) = length(levels.pressure)

"""$(TYPEDSIGNATURES)
Name of the vertical dimension in the output file or store for 3D atmospheric variables."""
vertical_dimension(::ModelLevels) = "layer"
vertical_dimension(::PressureLevels) = "pressure"

"""$(TYPEDSIGNATURES)
Vertical dimension of `variable` in `output`. Variables on the model's atmospheric layers
follow the output writer's `levels`, e.g. "layer" or "pressure". Variables with their own
vertical dimension (soil variables, custom output variables that extend the 1-argument
[`vertical_dimension`](@ref)) keep theirs."""
function vertical_dimension(output::AbstractOutput, variable::AbstractOutputVariable)
    dim = vertical_dimension(variable)
    return dim == "layer" ? vertical_dimension(output_levels(output)) : dim
end

"""$(TYPEDSIGNATURES)
Allocate the scratch field and the pressure levels on the model's architecture.
No-op for [`ModelLevels`](@ref)."""
initialize!(::ModelLevels, ::AbstractModel) = nothing

function initialize!(levels::PressureLevels, model::AbstractModel)
    # the interpolation itself doesn't care about the order (every level is searched for
    # independently) but a non-monotonic vertical coordinate in the output file would be
    # a mistake, so require monotonic, ascending (top to surface) or descending
    monotonic = issorted(levels.pressure) || issorted(levels.pressure, rev = true)
    @assert monotonic "Output pressure levels have to be monotonic, got $(levels.pressure)"
    @assert all(>(0), levels.pressure) "Output pressure levels have to be positive, got $(levels.pressure)"

    (; NF, grid) = model.spectral_grid
    arch = architecture(grid)

    # in the model's number format and on its architecture so that the interpolation
    # kernel can run on the GPU without scalar indexing into a host vector
    levels.pressure_on_arch = on_architecture(arch, NF.(levels.pressure))
    levels.scratch = Field(NF, grid, length(levels.pressure))
    return nothing
end

"""$(TYPEDSIGNATURES)
Vertical output levels of `output`, [`ModelLevels`](@ref) for output writers that don't
define a `levels` field."""
output_levels(output::AbstractOutput) = hasproperty(output, :levels) ? output.levels : ModelLevels()

"""$(TYPEDSIGNATURES)
Extrapolation of `variable` beyond the model's vertical levels: its `extrapolation` field
if it has one (e.g. [`TemperatureOutput`](@ref) descends dry-adiabatically below the
lowest model level), otherwise constant. Same `hasproperty` pattern as `transform` and
`unscale` in the generic [`output!`](@ref)."""
output_extrapolation(variable::AbstractOutputVariable) =
    hasproperty(variable, :extrapolation) ? variable.extrapolation : ConstantExtrapolation()

"""$(TYPEDSIGNATURES)
Interpolate `field` onto the output writer's vertical `levels`, returning the field to
write out. No-op on [`ModelLevels`](@ref) and for variables that are not 3D atmospheric
variables (2D variables, soil variables). For [`PressureLevels`](@ref) this interpolates
into the levels' scratch field, on the model grid and the model's architecture, i.e.
before the horizontal interpolation onto the output grid."""
interpolate_levels!(::ModelLevels, field, ::AbstractOutputVariable, ::AbstractSimulation) = field

function interpolate_levels!(
        levels::PressureLevels,
        field::AbstractField,
        variable::AbstractOutputVariable,
        simulation::AbstractSimulation,
    )
    # only 3D atmospheric variables, soil variables keep their own vertical dimension
    (is3D(variable) && !is_land(variable)) || return field

    (; model) = simulation
    surface_pressure = simulation.variables.parameterizations.surface_pressure

    interpolate_pressure_levels!(
        levels.scratch, field, surface_pressure, levels.pressure_on_arch,
        model.geometry.vertical_coordinates, levels.interpolation,
        output_extrapolation(variable),
    )
    return levels.scratch
end

"""$(TYPEDSIGNATURES)
Set κ of any [`DryAdiabaticExtrapolation`](@ref) in `output.variables` from the model's
atmosphere, so that the adiabatic descent below the lowest model level uses the same
κ = R_dry/cₚ as the model itself. Called at initialize!."""
function sync_extrapolations!(output::AbstractOutput, model::AbstractModel)
    hasproperty(model, :atmosphere) || return nothing
    for variable in values(output.variables)
        hasproperty(variable, :extrapolation) || continue
        variable.extrapolation = sync_extrapolation(variable.extrapolation, model.atmosphere)
    end
    return nothing
end

sync_extrapolation(E::AbstractVerticalExtrapolation, ::AbstractAtmosphere) = E
sync_extrapolation(::DryAdiabaticExtrapolation, atmosphere::AbstractAtmosphere) =
    DryAdiabaticExtrapolation(atmosphere.κ)
sync_extrapolation(E::SubsurfaceMask, atmosphere::AbstractAtmosphere) =
    SubsurfaceMask(sync_extrapolation(E.above_surface, atmosphere), E.missing_value)

"""$(TYPEDSIGNATURES)
Define the vertical coordinate of 3D atmospheric variables in the output file or store
`dest`: sigma for [`ModelLevels`](@ref), pressure in hPa for [`PressureLevels`](@ref).
Backend-agnostic via [`define_coordinate!`](@ref)."""
function define_vertical_coordinate!(dest, levels::ModelLevels, model::AbstractModel)
    σ = on_architecture(CPU(), model.geometry.σ_levels_full)
    return define_coordinate!(
        dest, vertical_dimension(levels), collect(Float64.(σ)),
        attribs = Dict("units" => "1", "long_name" => "sigma layer"),
    )
end

function define_vertical_coordinate!(dest, levels::PressureLevels, ::AbstractModel)
    return define_coordinate!(
        dest, vertical_dimension(levels), collect(Float64.(levels.pressure) ./ 100),
        attribs = Dict(
            "units" => "hPa", "long_name" => "pressure", "positive" => "down",
            "standard_name" => "air_pressure",
        ),
    )
end
