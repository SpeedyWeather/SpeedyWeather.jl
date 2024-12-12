abstract type AbstractTracer <: AbstractModelComponent end

export Tracer
@kwdef mutable struct Tracer <: AbstractTracer
    name::Symbol
    id::Int64 = 0
    active::Bool = true
end

Tracer(name::Symbol; kwargs...) = Tracer(; name, kwargs...)

const TRACER_DICT = Dict{Symbol, AbstractTracer}

function add!(dict::TRACER_DICT, tracers::Tracer...)
    for tracer in tracers
        if tracer.name in keys(dict)
            @warn "Tracer $(tracer.name) already exists. Skipping."
        else
            ntracers = length(dict)
            tracer.id = ntracers+1
            dict[tracer.name] = tracer
            @info "Tracer $(tracer.name) added as tracer $(tracer.id)"
        end
    end
end

add!(model::AbstractModel, tracers::Tracer...) = add!(model.tracers, tracers...)
