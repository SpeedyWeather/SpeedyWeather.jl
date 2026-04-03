# Setup function

"""Get nlayers for a given model type."""
nlayers_for_model(::Type{BarotropicModel}) = 1
nlayers_for_model(::Type{ShallowWaterModel}) = 1
nlayers_for_model(::Type{PrimitiveDryModel}) = 8
nlayers_for_model(::Type{PrimitiveWetModel}) = 8

"""Create a CPU model of the given type."""
function create_cpu_model(ModelType::Type; trunc = TRUNC, kwargs...)
    nlayers = nlayers_for_model(ModelType)
    spectral_grid = SpectralGrid(; nlayers, trunc)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = InitialConditions(spectral_grid, ModelType)
    return ModelType(spectral_grid; spectral_transform = M, feedback = nothing, initial_conditions, kwargs...)
end

"""Create a Reactant model of the given type."""
function create_reactant_model(ModelType::Type; trunc = TRUNC, kwargs...)
    nlayers = nlayers_for_model(ModelType)
    arch = SpeedyWeather.ReactantDevice()
    spectral_grid = SpectralGrid(; architecture = arch, nlayers, trunc)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = InitialConditions(spectral_grid, ModelType)
    return ModelType(spectral_grid; spectral_transform = M, feedback = nothing, initial_conditions, kwargs...)
end
