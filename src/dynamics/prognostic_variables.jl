function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

export PrognosticVariablesOcean
@kwdef struct PrognosticVariablesOcean{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridType,               # <: AbstractGrid
    GridVariable2D,         # <: AbstractField
} <: AbstractPrognosticVariables
    # DIMENSION
    "Grid used for ocean variables"
    grid::GridType

    # OCEAN VARIABLES
    "Sea surface temperature [K]"
    sea_surface_temperature::GridVariable2D = zeros(GridVariable2D, grid)

    "Sea ice concentration [1]"
    sea_ice_concentration::GridVariable2D = zeros(GridVariable2D, grid)

    "Prescribed ocean sensible heat flux [W/m²]"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Prescribed ocean evaporative flux [kg/s/m²]"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, grid)
end

export PrognosticVariablesLand
@kwdef struct PrognosticVariablesLand{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridType,               # <: AbstractGrid
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
} <: AbstractPrognosticVariables

    # DIMENSION
    "Grid used for land variables"
    grid::GridType

    "Number of soil layers for temperature and humidity"
    nlayers::Int

    # LAND VARIABLES
    "Soil temperature [K]"
    soil_temperature::GridVariable3D = zeros(GridVariable3D, grid, nlayers)

    "Soil moisture, volume fraction [1]"
    soil_moisture::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Snow depth [m]"
    snow_depth::GridVariable2D = zeros(GridVariable2D, grid)

    "Prescribed land sensible heat flux [W/m²], zero if not used"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Prescribed land evaporative flux [kg/s/m²], zero if not used"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, grid)
end

export PrognosticVariables
@kwdef struct PrognosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <: AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    SpectralVariable4D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    ParticleVector,         # <: AbstractVector{Particle{NF}}
    BoolVector,             # <: AbstractVector{Bool}   
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
    ocean::PrognosticVariablesOcean{NF, ArrayType, GridType, GridVariable2D} =
        PrognosticVariablesOcean{NF, ArrayType, GridType, GridVariable2D}(; grid)
    
    "Land variables, soil temperature, snow, and soil moisture"
    land::PrognosticVariablesLand{NF, ArrayType, GridType, GridVariable2D, GridVariable3D} =
        PrognosticVariablesLand{NF, ArrayType, GridType, GridVariable2D, GridVariable3D}(; grid, nlayers=nlayers_soil)

    "Tracers, last dimension is for n tracers [?]"
    tracers::Dict{Symbol, SpectralVariable4D} = Dict{Symbol, SpectralVariable4D}()

    "Particles for particle advection"
    particles::ParticleVector = zeros(ParticleVector, nparticles)

    "Particles activity"
    particles_activity::BoolVector = fill(true, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::Base.RefValue{NF} = Ref(one(NF))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::Clock = Clock()
end

Base.eltype(progn::PrognosticVariables{T}) where T = T

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
function PrognosticVariables(SG::SpectralGrid{Architecture, SpectrumType, GridType}; nsteps=DEFAULT_NSTEPS) where {Architecture, SpectrumType, GridType}

    (; spectrum, grid, nlayers, nlayers_soil, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector, BoolVector) = SG
    
    return PrognosticVariables{NF, ArrayType, SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector, BoolVector}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps,
        )
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid, model::AbstractModel)
    progn = PrognosticVariables(SG, nsteps = model.time_stepping.nsteps)
    add!(progn, model.tracers)
    return progn
end

function Base.show(
    io::IO,
    progn::PrognosticVariables{NF, ArrayType, SpectrumType, GridType},
) where {NF, ArrayType, SpectrumType, GridType}
    
    NFspectral = eltype(progn.vor)

    # resolution
    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn

    trunc = truncation(spectrum)
    Grid = RingGrids.nonparametric_type(GridType)
    nlat = RingGrids.get_nlat(grid)
    
    tracer_names = [key for (key, value) in progn.tracers]

    println(io, "PrognosticVariables{$NF, $ArrayType}")
    println(io, "├ vor:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ pres:  T$trunc, 1-layer, $nsteps-steps LowerTriangularArray{$NFspectral}")
    println(io, "├ random_pattern: T$trunc, 1-layer LowerTriangularArray{$NFspectral}")
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sea_ice_concentration:    Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       Field{$NF} on $nlat-ring $Grid")
    println(io, "│└ evaporative_flux:         Field{$NF} on $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
    println(io, "│├ soil_temperature:         Field{$NF} on $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ soil_moisture:            Field{$NF} on $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ snow_depth:               Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       Field{$NF} on $nlat-ring $Grid")
    println(io, "│└ evaporative_flux:         Field{$NF} on $nlat-ring $Grid")
    println(io, "├ tracers: $(length(tracer_names)), $tracer_names")
    println(io, "├ particles: $nparticles-element $(typeof(progn.particles))")
    println(io, "├ scale: $(progn.scale[])")
    print(io,   "└ clock: $(progn.clock.time)")
end

"""$(TYPEDSIGNATURES)
Copies entries of `progn_old` into `progn_new`."""
function Base.copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

    # Core variables using broadcast
    progn_new.vor .= progn_old.vor
    progn_new.div .= progn_old.div
    progn_new.temp .= progn_old.temp
    progn_new.humid .= progn_old.humid
    progn_new.pres .= progn_old.pres

    # Ocean variables
    progn_new.ocean.sea_surface_temperature .= progn_old.ocean.sea_surface_temperature
    progn_new.ocean.sea_ice_concentration .= progn_old.ocean.sea_ice_concentration
    progn_new.ocean.sensible_heat_flux .= progn_old.ocean.sensible_heat_flux
    progn_new.ocean.evaporative_flux .= progn_old.ocean.evaporative_flux

    # Land variables
    progn_new.land.soil_temperature .= progn_old.land.soil_temperature
    progn_new.land.snow_depth .= progn_old.land.snow_depth
    progn_new.land.soil_moisture .= progn_old.land.soil_moisture
    progn_new.land.sensible_heat_flux .= progn_old.land.sensible_heat_flux
    progn_new.land.evaporative_flux .= progn_old.land.evaporative_flux

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
        NF, ArrayType, SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector,
        }) where {
        NF, ArrayType, SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector,
        }

    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn
    
    # initialize regular progn variables 
    progn_new = PrognosticVariables{NF, ArrayType, SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, 
        ParticleVector}(;
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
    progn.ocean.evaporative_flux .= value

    # land
    progn.land.soil_temperature .= value
    progn.land.snow_depth .= value
    progn.land.soil_moisture .= value
    progn.land.sensible_heat_flux .= value
    progn.land.evaporative_flux .= value

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

"""$(TYPEDSIGNATURES)
Add `tracers` to the prognostic variables `progn` in `progn.tracers::Dict`."""
function add!(
    progn::PrognosticVariables{NF, ArrayType, SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D},
    tracers::Tracer...
) where {
        NF,                     # number format
        ArrayType,
        SpectrumType,
        GridType,
        SpectralVariable2D,
        SpectralVariable3D,
        SpectralVariable4D,
    }
    (; spectrum, nlayers, nsteps) = progn
    for tracer in tracers
        progn.tracers[tracer.name] = zeros(SpectralVariable4D, spectrum, nlayers, nsteps)
    end
end

"""$(TYPEDSIGNATURES)
Delete `tracers` in the prognostic variables `progn` in `progn.tracers::Dict`."""
function Base.delete!(progn::PrognosticVariables, tracers::Tracer...)
    for tracer in tracers
        delete!(progn.tracers, tracer.name)
    end
end