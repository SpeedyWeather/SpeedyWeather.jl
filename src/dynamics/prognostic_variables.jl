function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

# to be removed
struct PrognosticLayerTimesteps end
struct PrognosticSurfaceTimesteps end
struct PrognosticVariablesLayer end

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
end

export PrognosticVariablesLand
@kwdef struct PrognosticVariablesLand{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # LAND VARIABLES
    "Land surface temperature [K]"
    land_surface_temperature::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Snow depth [m]"
    snow_depth::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Soil moisture layer 1, volume fraction [1]"
    soil_moisture_layer1::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Soil moisture layer 2, volume fraction [1]"
    soil_moisture_layer2::GridVariable2D = zeros(GridVariable2D, nlat_half)
end

export PrognosticVariables
@kwdef struct PrognosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    NSTEPS,                 # number of timesteps
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    ParticleVector,         # <: AbstractVector{Particle{NF}}
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "max degree of spherical harmonics (0-based)"
    trunc::Int

    "Number of latitude rings on one hemisphere (Equator excl.), resolution parameter of grids"
    nlat_half::Int

    "number of vertical layers"
    nlayers::Int

    "Number of particles for particle advection"
    nparticles::Int

    # LAYERED VARIABLES
    "Vorticity of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    vor::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Absolute temperature [K]"
    temp::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Specific humidity [kg/kg]"
    humid::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Logarithm of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::NTuple{NSTEPS, SpectralVariable2D} =
        ntuple(i -> zeros(SpectralVariable2D, trunc+2, trunc+1), NSTEPS)

    "Ocean variables, sea surface temperature and sea ice concentration"
    ocean::PrognosticVariablesOcean{NF, ArrayType, GridVariable2D} =
        PrognosticVariablesOcean{NF, ArrayType, GridVariable2D}(; nlat_half)
    
    "Land variables, land surface temperature, snow and soil moisture"
    land::PrognosticVariablesLand{NF, ArrayType, GridVariable2D} =
        PrognosticVariablesLand{NF, ArrayType, GridVariable2D}(; nlat_half)

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
    (; trunc, nlat_half, nlayers, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector) = SG

    return PrognosticVariables{NF, ArrayType, nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector}(;
            trunc, nlat_half, nlayers=nlayers, nparticles,
        )
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid, model::AbstractModel)
    PrognosticVariables(SG, nsteps = model.time_stepping.nsteps)
end

function Base.show(
    io::IO,
    progn::PrognosticVariables{NF, ArrayType, NSTEPS},
) where {NF, ArrayType, NSTEPS}
    Grid = typeof(progn.ocean.sea_surface_temperature)
    println(io, "PrognosticVariables{$NF, $ArrayType}")
    
    # variables
    (; trunc, nlat_half, nlayers, nparticles) = progn
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    println(io, "├ vor:   T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ pres:  T$trunc, 1-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  $nlat-ring $Grid")
    println(io, "│└ sea_ice_concentration:    $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
    println(io, "│├ land_surface_temperature: $nlat-ring $Grid")
    println(io, "│├ soil_moisture_layer1:     $nlat-ring $Grid")
    println(io, "│└ soil_moisture_layer2:     $nlat-ring $Grid")
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
    progn_new.land.soil_moisture_layer1 .= progn_old.land.soil_moisture_layer1
    progn_new.land.soil_moisture_layer2 .= progn_old.land.soil_moisture_layer2

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