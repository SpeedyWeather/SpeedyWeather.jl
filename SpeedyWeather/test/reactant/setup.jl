# Setup function

"""Get nlayers for a given model type."""
nlayers_for_model(::Type{BarotropicModel}) = 1
nlayers_for_model(::Type{ShallowWaterModel}) = 1
nlayers_for_model(::Type{PrimitiveDryModel}) = 8
nlayers_for_model(::Type{PrimitiveWetModel}) = 8

"""Get default initial conditions for a given model type."""
function default_initial_conditions(::Type{BarotropicModel})
    return InitialConditions(; vordiv = ZeroInitially())
end

function default_initial_conditions(::Type{ShallowWaterModel})
    return InitialConditions(; vordiv = ZonalJet())
end

function default_initial_conditions(::Type{PrimitiveDryModel})
    return InitialConditions(; vordiv = ZonalWind())
end

function default_initial_conditions(::Type{PrimitiveWetModel})
    return InitialConditions(; vordiv = ZonalWind())
end

"""Create a CPU model of the given type."""
function create_cpu_model(ModelType::Type; trunc = TRUNC, kwargs...)
    nlayers = nlayers_for_model(ModelType)
    spectral_grid = SpectralGrid(; nlayers, trunc)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = default_initial_conditions(ModelType)
    return ModelType(spectral_grid; spectral_transform = M, initial_conditions, kwargs...)
end

"""Create a Reactant model of the given type."""
function create_reactant_model(ModelType::Type; trunc = TRUNC, kwargs...)
    nlayers = nlayers_for_model(ModelType)
    arch = SpeedyWeather.ReactantDevice()
    spectral_grid = SpectralGrid(; architecture = arch, nlayers, trunc)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = default_initial_conditions(ModelType)
    return ModelType(spectral_grid; spectral_transform = M, feedback = nothing, initial_conditions, kwargs...)
end
