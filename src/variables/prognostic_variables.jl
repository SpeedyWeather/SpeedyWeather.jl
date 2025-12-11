function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

export PrognosticVariables
@kwdef struct PrognosticVariables{
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <: AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    SpectralVariable4D,     # <: LowerTriangularArray
    OceanTuple,             # <: NamedTuple{Tuple{Symbol}, Tuple{PrognosticVariablesOcean}} #TODO: should the parameters change?
    LandTuple,              # <: NamedTuple{Tuple{Symbol}, Tuple{PrognosticVariablesLand}}
    PhysicsTuple,           # <: NamedTuple{Tuple{Symbol}, Tuple{PrognosticVariablesPhysics}}
    TracerTuple,            # <: NamedTuple{Tuple{Symbol}, Tuple{SpectralVariable4D}}
    ParticleVector,         # <: AbstractVector{Particle{NF}}
    RefValueNF,             # <: Base.RefValue{NF}
    ClockType,              # <: Union{Clock, Nothing}
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "spectral resolution"
    spectrum::SpectrumType

    "Grid used for variables"
    grid::GridType

    "number of vertical layers in the atmosphere"
    nlayers::Int

    "number of vertical layers in the soil"
    nlayers_soil::Int

    "Number of particles for particle advection"
    nparticles::Int

    "Number of time steps simultaneously stored in prognostic variables, 2 for 2-step leapfrog scheme"
    nsteps::Int

    # LAYERED VARIABLES
    "Vorticity of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    vor::SpectralVariable4D = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::SpectralVariable4D = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)

    "Absolute temperature [K]"
    temp::SpectralVariable4D = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)

    "Specific humidity [kg/kg]"
    humid::SpectralVariable4D = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)

    "Logarithm of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nsteps)

    "Random pattern following a random process [1]"
    random_pattern::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)

    "Ocean variables, sea surface temperature and sea ice concentration"
    ocean::OceanTuple = NamedTuple()
    
    "Land variables, soil temperature, snow, and soil moisture"
    land::LandTuple = NamedTuple()

    "Physics variables from the parametrizations"
    physics::PhysicsTuple = NamedTuple()
    
    "Tracers, last dimension is for n tracers [?]"
    tracers::TracerTuple = NamedTuple()

    "Particles for particle advection"
    particles::ParticleVector = zeros(ParticleVector, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::RefValueNF = Ref(one(real(eltype(vor))))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::ClockType = Clock()
end

Adapt.@adapt_structure PrognosticVariables

Base.eltype(::PrognosticVariables{SpectrumType, GridType, SpectralVariable2D}) where {SpectrumType, GridType, SpectralVariable2D <: LowerTriangularArray{Complex{NF}}} where NF = NF
Architectures.array_type(::PrognosticVariables{SpectrumType, GridType, SpectralVariable2D}) where {SpectrumType, GridType, SpectralVariable2D <: LowerTriangularArray{NF, N, ArrayType}} where {NF, N, ArrayType <: AbstractArray} = nonparametric_type(ArrayType)

function get_steps(coeffs::LowerTriangularArray{T, 2}) where T
    nsteps = size(coeffs, 2)
    return ntuple(i -> lta_view(coeffs, :, i), nsteps)
end

function get_steps(coeffs::LowerTriangularArray{T, 3}) where T
    nsteps = size(coeffs, 3)
    return ntuple(i -> lta_view(coeffs, :, :, i), nsteps)
end

export get_step

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables used for the leapfrog time step.
This method is for a 2D spectral variable (horizontal only) with steps in the 3rd dimension."""
get_step(coeffs::LowerTriangularArray{T, 2}, i) where T = lta_view(coeffs, :, i)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables used for the leapfrog time step.
This method is for a 3D spectral variable (horizontal+vertical) with steps in the 4rd dimension."""
get_step(coeffs::LowerTriangularArray{T, 3}, i) where T = lta_view(coeffs, :, :, i)

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(model::AbstractModel)

    SG = model.spectral_grid 
    tracers = model.tracers
    nsteps = model.time_stepping.nsteps

    (; NF, spectrum, grid, nlayers, nlayers_soil, nparticles) = SG
    (; SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, ParticleVector) = SG
    
    # allocate parameterization variables 
    variable_names = get_prognostic_variables(model)

    # TODO: currently this is just a drop-in replacement, later we should have nlayers_soil from the land model and not from the SpectralGrid
    land = initialize_variables(SG, nlayers_soil, variable_names.land...)
    physics = initialize_variables(SG, nlayers, variable_names.atmosphere...)
    ocean = initialize_variables(SG, 1, variable_names.ocean...)

    tracer_tuple = (; [key => zeros(SpectralVariable4D, spectrum, nlayers, nsteps) for key in keys(tracers)]...)
    clock = Clock()

    return PrognosticVariables{typeof(spectrum), typeof(grid),
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, typeof(ocean), typeof(land), typeof(physics), typeof(tracer_tuple), ParticleVector, Base.RefValue{NF}, typeof(clock)}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps, ocean, land, physics, tracers = tracer_tuple, clock
        )
end

function Base.show(
    io::IO,
    progn::PrognosticVariables{SpectrumType, GridType},
) where {SpectrumType, GridType}
    
    NF = eltype(progn)
    NFspectral = eltype(progn.vor)
    ArrayType = array_type(progn)

    # resolution
    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn

    trunc = truncation(spectrum)
    Grid = nonparametric_type(GridType)
    nlat = RingGrids.get_nlat(grid)
    
    tracer_names = keys(progn.tracers)

    println(io, "PrognosticVariables{$NF, $ArrayType}")
    println(io, "├ vor:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ pres:  T$trunc, 1-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ random_pattern: T$trunc, 1-layer LowerTriangularArray{$NFspectral}")
    println(io, "├┐ocean:")
    println(io, "│├ sea_surface_temperature:  Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sea_ice_concentration:    Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       Field{$NF} on $nlat-ring $Grid")
    println(io, "│└ surface_humidity_flux:    Field{$NF} on $nlat-ring $Grid")
    println(io, "├┐land:")
    println(io, "│├ soil_temperature:         Field{$NF} on $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ soil_moisture:            Field{$NF} on $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ snow_depth:               Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       Field{$NF} on $nlat-ring $Grid")
    println(io, "│└ surface_humidity_flux:    Field{$NF} on $nlat-ring $Grid")
    println(io, "├ tracers: $(length(tracer_names)), $tracer_names")
    println(io, "├ particles: $nparticles-element $(typeof(progn.particles))")
    println(io, "├ scale: $(progn.scale[])")
    print(io,   "└ clock: $(progn.clock.time)")
end

function copy_if_key_exists!(to, from, key)
    # use hasproperty here as a union for haskey (works with NamedTuples) and hasfield (works with structs)
    if hasproperty(to, key) && hasproperty(from, key)
        getfield(to, key) .= getfield(from, key)
    end
end

"""$(TYPEDSIGNATURES)
Copies entries of `progn_old` into `progn_new`."""
function Base.copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

    # Atmospheric variables
    copy_if_key_exists!(progn_new, progn_old, :vor)
    copy_if_key_exists!(progn_new, progn_old, :div)
    copy_if_key_exists!(progn_new, progn_old, :temp)
    copy_if_key_exists!(progn_new, progn_old, :humid)
    copy_if_key_exists!(progn_new, progn_old, :pres)

    # Ocean variables
    copy_if_key_exists!(progn_new.ocean, progn_old.ocean, :sea_surface_temperature)
    copy_if_key_exists!(progn_new.ocean, progn_old.ocean, :sea_ice_concentration)
    copy_if_key_exists!(progn_new.ocean, progn_old.ocean, :sensible_heat_flux)
    copy_if_key_exists!(progn_new.ocean, progn_old.ocean, :surface_humidity_flux)

    # Land variables
    copy_if_key_exists!(progn_new.land, progn_old.land, :soil_temperature)
    copy_if_key_exists!(progn_new.land, progn_old.land, :snow_depth)
    copy_if_key_exists!(progn_new.land, progn_old.land, :soil_moisture)
    copy_if_key_exists!(progn_new.land, progn_old.land, :sensible_heat_flux)
    copy_if_key_exists!(progn_new.land, progn_old.land, :surface_humidity_flux)

    # Tracers - using broadcast assignment
    for (key, value) in progn_old.tracers
        progn_new.tracers[key] .= value
    end

    # Random pattern
    progn_new.random_pattern .= progn_old.random_pattern

    # Particles - copy only up to the minimum length
    n = min(length(progn_new.particles), length(progn_old.particles))
    if n > 0
        progn_new.particles[1:n] .= progn_old.particles[1:n]
    end

    # Copy scale and clock
    copy!(progn_new.clock, progn_old.clock)
    progn_new.scale[] = progn_old.scale[]

    return nothing
end

function Base.zero(
    progn::PrognosticVariables{
        SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, RefValueNF, ClockType,
        }) where {
        SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, RefValueNF, ClockType,
        }

    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn
    
    # initialize regular progn variables 
    progn_new = PrognosticVariables{SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple,
        ParticleVector, RefValueNF, ClockType}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps
        )

    # add tracers with zero 
    for (key, value) in progn.tracers 
        progn_new.tracers[key] = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)
    end 

    # use the same scale 
    progn_new.scale[] = progn.scale[]

    return progn_new
end 

function Base.fill!(progn::PrognosticVariables, value::Number)

    progn.vor .= value
    progn.div .= value
    progn.temp .= value
    progn.humid .= value
    progn.pres .= value

    #TODO copy over random pattern?

    # ocean
    progn.ocean.sea_surface_temperature .= value
    progn.ocean.sea_ice_concentration .= value
    progn.ocean.sensible_heat_flux .= value
    progn.ocean.surface_humidity_flux .= value

    # land
    progn.land.soil_temperature .= value
    progn.land.snow_depth .= value
    progn.land.soil_moisture .= value
    progn.land.sensible_heat_flux .= value
    progn.land.surface_humidity_flux .= value

    # fill tracers
    for (key, value) in progn.tracers 
        for value_i in value # istep of nsteps tuple 
            value_i .= value
        end 
    end 

    # particles are ignored for the fill

    return progn
end 

function Base.one(progn::PrognosticVariables)
    zero_progn = zero(progn)
    fill!(zero_progn, 1)
    return zero_progn
end 
