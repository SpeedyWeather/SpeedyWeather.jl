export WriteVariablesRestartFile

"""WriteVariablesRestartFile callback. Writes a restart file for the prognostic variables
on the last time step as jld2 file. Options are $(TYPEDFIELDS)"""
@kwdef mutable struct WriteVariablesRestartFile <: AbstractCallback
    "[OPTION] File name for restart file, should end with .jld2"
    filename::String = "restart.jld2"

    "[OPTION] Path for restart file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] Apply lossless compression in JLD2?"
    compress::Bool = true

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true

    "[DERIVED] package version to track possible incompatibilities"
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)
end

initialize!(::WriteVariablesRestartFile, args...) = nothing
callback!(::WriteVariablesRestartFile, args...) = nothing

"""
$(TYPEDSIGNATURES)
A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions.
Variables in restart file are not scaled."""
function finalize!(
        restart::WriteVariablesRestartFile,
        vars::Variables,
        model::AbstractModel,
    )
    # escape in case of no output
    restart.write_only_with_output && (model.output.active || return nothing)

    (; compress, filename) = restart

    # use output run path if not specified
    path = restart.path == "" ? model.output.run_path : restart.path
    mkpath(path)
    restart.path = path     # update path in case it was empty before

    jldopen(joinpath(path, filename), "w"; compress) do f
        f["variables.prognostic"] = vars.prognostic
        f["version"] = restart.pkg_version
        f["description"] = "Restart file created by SpeedyWeather.jl"
    end

    return model.output.active || @info "Restart file written to $(joinpath(path, filename)) although output=false"
end

export WriteModelComponentFile

"""WriteModelComponentFile callback. Writes a file containing a model component as .jld2.
This is written at the end of the simulation, reflecting therefore the state that model
component at the end of the simulation. Options are $(TYPEDFIELDS)"""
@kwdef mutable struct WriteModelComponentFile{C} <: AbstractCallback
    "[OPTION] File name for model component restart file, should end with .jld2"
    filename::String = "model_component.jld2"

    "[OPTION] Path for restart file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] Apply lossless compression in JLD2?"
    compress::Bool = true

    "[OPTION] Only write with model.output.active = true?"
    write_only_with_output::Bool = true

    "[OPTION] Model component to write, to be passed on as keyword argument to the callback"
    component::C

    "[DERIVED] package version to track possible incompatibilities"
    pkg_version::VersionNumber = isnothing(pkgversion(SpeedyWeather)) ? v"0.0.0" : pkgversion(SpeedyWeather)
end

initialize!(::WriteModelComponentFile, args...) = nothing
callback!(::WriteModelComponentFile, args...) = nothing

"""$(TYPEDSIGNATURES)
A restart file for the specified model component is written
to the output folder (or current path) that can be used to restart the model.
With the `component` field in `WriteModelComponentFile` as at the end of the simulation."""
function finalize!(
        restart::WriteModelComponentFile,
        ::Variables,
        model::AbstractModel,
    )
    # escape in case of no output
    restart.write_only_with_output && (model.output.active || return nothing)

    (; compress, filename) = restart

    # use output run path if not specified
    path = restart.path == "" ? model.output.run_path : restart.path
    mkpath(path)
    restart.path = path     # update path in case it was empty before

    jldopen(joinpath(path, filename), "w"; compress) do f
        f["component"] = restart.component
        f["version"] = restart.pkg_version
        f["description"] = "Model component file created by SpeedyWeather.jl"
    end

    return model.output.active || @info "Model component file written to $(joinpath(path, filename)) although output=false"
end

load_model_component(callback::WriteModelComponentFile) =
    load_model_component(callback.path, callback.filename)

function load_model_component(path::String, filename::String = "")
    jldopen(joinpath(path, filename), "r") do f
        return f["component"]
    end
end