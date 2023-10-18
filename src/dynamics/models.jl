"""
$(TYPEDSIGNATURES)
Simulation is a container struct to be used with `run!(::Simulation)`.
It contains
$(TYPEDFIELDS)"""
struct Simulation{Model<:ModelSetup}
    "define the current state of the model"
    prognostic_variables::PrognosticVariables

    "contain the tendencies and auxiliary arrays to compute them"
    diagnostic_variables::DiagnosticVariables

    "all parameters, constant at runtime"
    model::Model
end

"""
$(SIGNATURES)
The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: Barotropic
    "dictates resolution for many other components"
    spectral_grid::SpectralGrid = SpectralGrid(nlev=1)

    # DYNAMICS
    "contains physical and orbital characteristics"
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    forcing::AbstractForcing{NF} = NoForcing(spectral_grid)
    drag::AbstractDrag{NF} = NoDrag(spectral_grid)
    initial_conditions::InitialConditions = StartWithRandomVorticity()

    # NUMERICS
    time_stepping::TimeStepper{NF} = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit{NF} = NoImplicit(spectral_grid)

    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,Barotropic)
    feedback::AbstractFeedback = Feedback()
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!`."""
function initialize!(model::Barotropic)
    (;spectral_grid) = model

    spectral_grid.nlev > 1 && @warn "Only nlev=1 supported for BarotropicModel, \
        SpectralGrid with nlev=$(spectral_grid.nlev) provided."

    initialize!(model.forcing,model)
    initialize!(model.drag,model)
    initialize!(model.horizontal_diffusion,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end

"""
$(SIGNATURES)
The ShallowWaterModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ShallowWater
    "dictates resolution for many other components"
    spectral_grid::SpectralGrid = SpectralGrid(nlev=1)

    # DYNAMICS
    "contains physical and orbital characteristics"
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    forcing::AbstractForcing{NF} = NoForcing(spectral_grid)
    drag::AbstractDrag{NF} = NoDrag(spectral_grid)
    initial_conditions::InitialConditions = ZonalJet()
    orography::AbstractOrography{NF} = EarthOrography(spectral_grid)

    # NUMERICS
    time_stepping::TimeStepper{NF} = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit{NF} = ImplicitShallowWater(spectral_grid)

    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,ShallowWater)
    feedback::AbstractFeedback = Feedback()
end

has(::Type{<:ShallowWater}, var_name::Symbol) = var_name in (:vor, :div, :pres)
default_concrete_model(::Type{ShallowWater}) = ShallowWaterModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::ShallowWater)
    (;spectral_grid) = model

    spectral_grid.nlev > 1 && @warn "Only nlev=1 supported for ShallowWaterModel, \
                                SpectralGrid with nlev=$(spectral_grid.nlev) provided."

    initialize!(model.forcing,model)
    initialize!(model.drag,model)
    initialize!(model.horizontal_diffusion,model)
    initialize!(model.orography,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end

"""
$(SIGNATURES)
The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef struct PrimitiveDryModel{NF<:AbstractFloat, D<:AbstractDevice} <: PrimitiveDry
    "dictates resolution for many other components"
    spectral_grid::SpectralGrid = SpectralGrid()

    # DYNAMICS
    "contains physical and orbital characteristics"
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    initial_conditions::InitialConditions = ZonalWind()
    orography::AbstractOrography{NF} = EarthOrography(spectral_grid)

    # BOUNDARY CONDITIONS
    land_sea_mask::AbstractLandSeaMask{NF} = LandSeaMask(spectral_grid)
    ocean::AbstractOcean{NF} = SeasonalOceanClimatology(spectral_grid)
    land::AbstractLand{NF} = SeasonalLandClimatology(spectral_grid)

    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    boundary_layer_drag::BoundaryLayerDrag{NF} = LinearDrag(spectral_grid)
    temperature_relaxation::TemperatureRelaxation{NF} = HeldSuarez(spectral_grid)
    static_energy_diffusion::VerticalDiffusion{NF} = StaticEnergyDiffusion(spectral_grid)
    surface_thermodynamics::AbstractSurfaceThermodynamics{NF} = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::AbstractSurfaceWind{NF} = SurfaceWind(spectral_grid)
    surface_heat_flux::AbstractSurfaceHeat{NF} = SurfaceSensibleHeat(spectral_grid)
    
    # NUMERICS
    time_stepping::TimeStepper{NF} = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit{NF} = ImplicitPrimitiveEq(spectral_grid)
    vertical_advection::VerticalAdvection{NF} = CenteredVerticalAdvection(spectral_grid)
    
    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,PrimitiveDry)
    feedback::AbstractFeedback = Feedback()
end

has(::Type{<:PrimitiveDry}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
default_concrete_model(::Type{PrimitiveDry}) = PrimitiveDryModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveDry)
    (;spectral_grid) = model

    # numerics (implicit is initialized later)
    initialize!(model.horizontal_diffusion,model)

    # boundary conditionss
    initialize!(model.orography,model)
    initialize!(model.land_sea_mask)
    initialize!(model.ocean)
    initialize!(model.land)

    # parameterizations
    initialize!(model.boundary_layer_drag,model)
    initialize!(model.temperature_relaxation,model)
    initialize!(model.static_energy_diffusion,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)
    
    (;time) = prognostic_variables.clock
    initialize!(prognostic_variables.ocean,time,model)
    initialize!(prognostic_variables.land,time,model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end

"""
$(SIGNATURES)
The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef struct PrimitiveWetModel{NF<:AbstractFloat, D<:AbstractDevice} <: PrimitiveWet
    "dictates resolution for many other components"
    spectral_grid::SpectralGrid = SpectralGrid()

    # DYNAMICS
    "contains physical and orbital characteristics"
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    initial_conditions::InitialConditions = ZonalWind()
    orography::AbstractOrography{NF} = EarthOrography(spectral_grid)

    # BOUNDARY CONDITIONS
    land_sea_mask::AbstractLandSeaMask{NF} = LandSeaMask(spectral_grid)
    ocean::AbstractOcean{NF} = SeasonalOceanClimatology(spectral_grid)
    land::AbstractLand{NF} = SeasonalLandClimatology(spectral_grid)

    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    thermodynamics::Thermodynamics{NF} = Thermodynamics(spectral_grid,atmosphere)
    boundary_layer_drag::BoundaryLayerDrag{NF} = LinearDrag(spectral_grid)
    temperature_relaxation::TemperatureRelaxation{NF} = HeldSuarez(spectral_grid)
    static_energy_diffusion::VerticalDiffusion{NF} = StaticEnergyDiffusion(spectral_grid)
    large_scale_condensation::AbstractCondensation{NF} = SpeedyCondensation(spectral_grid)
    surface_thermodynamics::AbstractSurfaceThermodynamics{NF} = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::AbstractSurfaceWind{NF} = SurfaceWind(spectral_grid)
    surface_heat_flux::AbstractSurfaceHeat{NF} = SurfaceSensibleHeat(spectral_grid)
    evaporation::AbstractEvaporation{NF} = SurfaceEvaporation(spectral_grid)
    
    # NUMERICS
    time_stepping::TimeStepper{NF} = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit{NF} = ImplicitPrimitiveEq(spectral_grid)
    vertical_advection::VerticalAdvection{NF} = CenteredVerticalAdvection(spectral_grid)
    
    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,PrimitiveWet)
    feedback::AbstractFeedback = Feedback()
end
 
has(::Type{<:PrimitiveWet}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres, :humid)
default_concrete_model(::Type{PrimitiveWet}) = PrimitiveWetModel
 
"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveWet)
    (;spectral_grid) = model

    # numerics (implicit is initialized later)
    initialize!(model.horizontal_diffusion,model)

    # boundary conditionss
    initialize!(model.orography,model)
    initialize!(model.land_sea_mask)
    initialize!(model.ocean)
    initialize!(model.land)

    # parameterizations
    initialize!(model.boundary_layer_drag,model)
    initialize!(model.temperature_relaxation,model)
    initialize!(model.static_energy_diffusion,model)
    initialize!(model.large_scale_condensation,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)

    (;time) = prognostic_variables.clock
    initialize!(prognostic_variables.ocean,time,model)
    initialize!(prognostic_variables.land,time,model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end
   
"""$(TYPEDSIGNATURES)
Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
The default fallback is that all variables are included. """
has(M::Type{<:ModelSetup}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
has(M::ModelSetup, var_name) = has(typeof(M), var_name)

# strip away the parameters of the model type
model_class(::Type{<:Barotropic}) = Barotropic
model_class(::Type{<:ShallowWater}) = ShallowWater
model_class(::Type{<:PrimitiveDry}) = PrimitiveDry
model_class(::Type{<:PrimitiveWet}) = PrimitiveWet
model_class(model::ModelSetup) = model_class(typeof(model))

function Base.show(io::IO,M::ModelSetup)
    println(io,"$(typeof(M))")
    for key in propertynames(M)[1:end-1]
        val = getfield(M,key)
        println(io,"├ $key: $(typeof(val))")
    end
    print(io,"└ feedback: $(typeof(M.feedback))")
end

function Base.show(io::IO,S::Simulation)
    println(io,"$(typeof(S))")
    println(io,"├ $(typeof(S.model))")
    println(io,"├ $(typeof(S.prognostic_variables))")
    print(io,  "└ $(typeof(S.diagnostic_variables))")
end