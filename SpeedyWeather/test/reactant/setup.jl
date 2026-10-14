# Setup functions

"""Get nlayers for a given model type."""
nlayers_for_model(::Type{BarotropicModel}) = 1
nlayers_for_model(::Type{ShallowWaterModel}) = 1
nlayers_for_model(::Type{PrimitiveDryModel}) = 8
nlayers_for_model(::Type{PrimitiveWetModel}) = 8

"""Model-specific keyword arguments for the model constructor."""
model_kwargs(::Type{<:Union{BarotropicModel, ShallowWaterModel}}, spectral_grid) = (;)

# PrimitiveDryModel already defaults to a constant longwave transmissivity
model_kwargs(::Type{PrimitiveDryModel}, spectral_grid) = (; convection = nothing)

function model_kwargs(::Type{PrimitiveWetModel}, spectral_grid)
    # constant transmissivity: Frierson's surface-concentrated absorption without convection
    # amplifies ULP-level CPU-vs-Reactant differences exponentially (chaotic, not a Reactant bug)
    transmissivity = ConstantLongwaveTransmissivity(spectral_grid)
    radiation = Radiation(spectral_grid, longwave = OneBandLongwave(spectral_grid; transmissivity))
    return (; radiation, convection = nothing)
end

"""Create a model of the given type on architecture `arch`."""
function create_model(ModelType::Type, arch; truncation = TRUNCATION, kwargs...)
    nlayers = nlayers_for_model(ModelType)
    spectral_grid = SpectralGrid(; architecture = arch, nlayers, truncation)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = InitialConditions(spectral_grid, ModelType)
    return ModelType(
        spectral_grid; spectral_transform = M, feedback = nothing, initial_conditions,
        model_kwargs(ModelType, spectral_grid)..., kwargs...
    )
end

"""Create a CPU model of the given type."""
create_cpu_model(ModelType::Type; kwargs...) = create_model(ModelType, SpeedyWeather.CPU(); kwargs...)

"""Create a Reactant model of the given type."""
create_reactant_model(ModelType::Type; kwargs...) = create_model(ModelType, SpeedyWeather.ReactantDevice(); kwargs...)

"""Synchronize prognostic variables from Reactant to CPU simulation
using copy_variables! which handles arrays (via copyto!), Ref values,
and mutable structs (Clock) recursively."""
function sync_variables!(sim_cpu, sim_reactant)
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_reactant, _ = SpeedyWeather.unpack(sim_reactant)

    # Copy all variables
    SpeedyWeather.copy!(vars_cpu, vars_reactant)
    return
end
