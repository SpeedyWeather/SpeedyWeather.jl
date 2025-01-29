function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

export PrognosticVariablesOcean
@kwdef struct PrognosticVariablesOcean{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # OCEAN VARIABLES
    "Sea surface temperature [K]"
    sea_surface_temperature::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Sea ice concentration [1]"
    sea_ice_concentration::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Prescribed ocean sensible heat flux [W/m²]"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Prescribed ocean evaporative flux [kg/s/m²]"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)
end

export PrognosticVariablesLand
@kwdef struct PrognosticVariablesLand{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,
} <: AbstractPrognosticVariables

    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    "Number of soil layers for temperature and humidity"
    nlayers::Int

    # LAND VARIABLES
    "Soil temperature [K]"
    soil_temperature::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)

    "Soil moisture, volume fraction [1]"
    soil_moisture::GridVariable3D = zeros(GridVariable2D, nlat_half, nlayers)
    
    "Snow depth [m]"
    snow_depth::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Prescribed land sensible heat flux [W/m²], zero if not used"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Prescribed land evaporative flux [kg/s/m²], zero if not used"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)
end

export PrognosticVariables
@kwdef struct PrognosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    nsteps,                 # number of timesteps
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
    ParticleVector,         # <: AbstractVector{Particle{NF}}
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "max degree of spherical harmonics (0-based)"
    trunc::Int

    "Number of latitude rings on one hemisphere (Equator excl.), resolution parameter of grids"
    nlat_half::Int

    "number of vertical layers in the atmosphere"
    nlayers::Int

    "number of vertical layers in the soil"
    nlayers_soil::Int

    "Number of particles for particle advection"
    nparticles::Int

    # LAYERED VARIABLES
    "Vorticity of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    vor::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), nsteps)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), nsteps)

    "Absolute temperature [K]"
    temp::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), nsteps)

    "Specific humidity [kg/kg]"
    humid::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), nsteps)

    "Logarithm of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::NTuple{nsteps, SpectralVariable2D} =
        ntuple(i -> zeros(SpectralVariable2D, trunc+2, trunc+1), nsteps)

    "Random pattern following a random process [1]"
    random_pattern::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)

    "Ocean variables, sea surface temperature and sea ice concentration"
    ocean::PrognosticVariablesOcean{NF, ArrayType, GridVariable2D} =
        PrognosticVariablesOcean{NF, ArrayType, GridVariable2D}(; nlat_half)
    
    "Land variables, soil temperature, snow, and soil moisture"
    land::PrognosticVariablesLand{NF, ArrayType, GridVariable2D, GridVariable3D} =
        PrognosticVariablesLand{NF, ArrayType, GridVariable2D, GridVariable3D}(; nlat_half, nlayers=nlayers_soil)

    "Tracers, last dimension is for n tracers [?]"
    tracers::Dict{Symbol, NTuple{nsteps, SpectralVariable3D}} = Dict{Symbol, NTuple{nsteps, SpectralVariable3D}}()

    "Particles for particle advection"
    particles::ParticleVector = zeros(ParticleVector, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::Base.RefValue{NF} = Ref(one(NF))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::Clock = Clock()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid; nsteps=DEFAULT_NSTEPS)

    (; trunc, nlat_half, nlayers, nlayers_soil, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector) = SG

    return PrognosticVariables{NF, ArrayType, nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector}(;
            trunc, nlat_half, nlayers, nlayers_soil, nparticles,
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
    progn::PrognosticVariables{NF, ArrayType, nsteps},
) where {NF, ArrayType, nsteps}
    
    Grid = typeof(progn.ocean.sea_surface_temperature)
    tracer_names = [key for (key, value) in progn.tracers]
    
    println(io, "PrognosticVariables{$NF, $ArrayType}")
    
    # variables
    (; trunc, nlat_half, nlayers, nlayers_soil, nparticles) = progn
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    println(io, "├ vor:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ pres:  T$trunc, 1-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ random_pattern: T$trunc, 1-layer LowerTriangularArray{$NF}")
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  $nlat-ring $Grid")
    println(io, "│└ sea_ice_concentration:    $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
    println(io, "│├ soil_temperature: $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ soil_moisture:    $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│└ snow_depth: $nlat-ring $Grid")
    println(io, "├ tracers: $(length(tracer_names)), $tracer_names")
    println(io, "├ particles: $nparticles-element $(typeof(progn.particles))")
    println(io, "├ scale: $(progn.scale[])")
    print(io,   "└ clock: $(progn.clock.time)")
end

# has(::PrognosticVariables{NF, Grid, M}, var_name::Symbol) where {NF, Grid, M} = has(M, var_name)

"""$(TYPEDSIGNATURES)
Copies entries of `progn_old` into `progn_new`."""
function Base.copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

    for i in eachindex(progn_new.vor)   # each leapfrog time step
        progn_new.vor[i] .= progn_old.vor[i]
        progn_new.div[i] .= progn_old.div[i]
        progn_new.temp[i] .= progn_old.temp[i]
        progn_new.humid[i] .= progn_old.humid[i]
        progn_new.pres[i] .= progn_old.pres[i]
    end

    # ocean
    progn_new.ocean.sea_surface_temperature .= progn_old.ocean.sea_surface_temperature
    progn_new.ocean.sea_ice_concentration .= progn_old.ocean.sea_ice_concentration
    
    # land
    progn_new.land.land_surface_temperature .= progn_old.land.land_surface_temperature
    progn_new.land.snow_depth .= progn_old.land.snow_depth
    progn_new.land.soil_moisture_layer1 .= progn_old.land.soil_moisture_layer1
    progn_new.land.soil_moisture_layer2 .= progn_old.land.soil_moisture_layer2

    # TODO copy over tracers

    # copy largest subset of particles
    if length(progn_new.particles) != length(progn_old.particles)
        nnew = length(progn_new.particles)
        nold = length(progn_old.particles)
        nsub = min(nnew, nold)
        @warn "Number of particles changed (origin: $nold, destination: $nnew), copying over only the largest subset ($nsub particles)"
        progn_new.particles[1:nsub] .= progn_old.particles[1:nsub]
    else
        progn_new.particles .= progn_old.particles
    end

    progn_new.clock.time = progn_old.clock.time
    progn_new.scale[] = progn_old.scale[]

    return progn_new
end

"""$(TYPEDSIGNATURES)
Add `tracers` to the prognostic variables `progn` in `progn.tracers::Dict`."""
function add!(progn::PrognosticVariables{NF, ArrayType, nsteps, SpectralVariable2D, SpectralVariable3D}, tracers::Tracer...
    ) where {
        NF,                     # number format
        ArrayType,
        nsteps,
        SpectralVariable2D,
        SpectralVariable3D,
    }
    (; trunc, nlayers) = progn
    for tracer in tracers
        progn.tracers[tracer.name] = ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), nsteps)
    end
end

"""$(TYPEDSIGNATURES)
Delete `tracers` in the prognostic variables `progn` in `progn.tracers::Dict`."""
function Base.delete!(progn::PrognosticVariables, tracers::Tracer...)
    for tracer in tracers
        delete!(progn.tracers, tracer.name)
    end
end