function Base.show(io::IO,C::AbstractCallback)
    println(io,"$(typeof(C)) <: AbstractCallback")
    keys = propertynames(C)
    print_fields(io,C,keys)
end

# dummy callback
export NoCallback
struct NoCallback <: AbstractCallback end
initialize!(::NoCallback,args...) = nothing
callback!(::NoCallback,args...) = nothing
finish!(::NoCallback,args...) = nothing

# simply loop over vector of callbacks
for func in (:initialize!, :callback!, :finish!)
    @eval begin
        function $func(callbacks::Vector{<:AbstractCallback},args...)
            for callback in callbacks
                $func(callback,args...)
            end
        end
    end
end

export GlobalSurfaceTemperatureCallback
Base.@kwdef mutable struct GlobalSurfaceTemperatureCallback{NF} <: AbstractCallback
    timestep_counter::Int = 0
    temp::Vector{NF} = [0]
end

GlobalSurfaceTemperatureCallback(SG::SpectralGrid) = GlobalSurfaceTemperatureCallback{SG.NF}()

function initialize!(
    callback::GlobalSurfaceTemperatureCallback{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
) where NF
    callback.temp = Vector{NF}(undef,progn.clock.n_timesteps+1)
    callback.temp[1] = diagn.layers[diagn.nlev].temp_average[]
    callback.timestep_counter = 1
end

function callback!(
    callback::GlobalSurfaceTemperatureCallback,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::ModelSetup,
)
    callback.timestep_counter += 1
    i = callback.timestep_counter
    callback.temp[i] = diagn.layers[diagn.nlev].temp_average[]
end

finish!(::GlobalSurfaceTemperatureCallback,args...) = nothing