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
    soil_moisture::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
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
    SpectrumType,           # <: AbstractSpectrum
    nsteps,                 # number of timesteps
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
    ParticleVector,         # <: AbstractVector{Particle{NF}}
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "spectral resolution"
    spectrum::SpectrumType

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
        ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)

    "Absolute temperature [K]"
    temp::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)

    "Specific humidity [kg/kg]"
    humid::NTuple{nsteps, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)

    "Logarithm of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::NTuple{nsteps, SpectralVariable2D} =
        ntuple(i -> zeros(SpectralVariable2D, spectrum), nsteps)

    "Random pattern following a random process [1]"
    random_pattern::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)

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

    (; spectrum, nlat_half, nlayers, nlayers_soil, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector) = SG

    return PrognosticVariables{NF, ArrayType, typeof(spectrum), nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector}(;
            spectrum, nlat_half, nlayers, nlayers_soil, nparticles,
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
    progn::PrognosticVariables{NF, ArrayType, SpectrumType, nsteps},
) where {NF, ArrayType, SpectrumType, nsteps}
    
    Grid = typeof(progn.ocean.sea_surface_temperature)
    tracer_names = [key for (key, value) in progn.tracers]
    
    println(io, "PrognosticVariables{$NF, $ArrayType}")
    
    # variables
    (; spectrum, nlat_half, nlayers, nlayers_soil, nparticles) = progn
    trunc = truncation(spectrum)
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    println(io, "├ vor:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ pres:  T$trunc, 1-layer, $nsteps-steps LowerTriangularArray{$NF}")
    println(io, "├ random_pattern: T$trunc, 1-layer LowerTriangularArray{$NF}")
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  $nlat-ring $Grid")
    println(io, "│├ sea_ice_concentration:    $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       $nlat-ring $Grid")
    println(io, "│└ evaporative_flux:         $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
    println(io, "│├ soil_temperature:         $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ soil_moisture:            $nlayers_soil-layer, $nlat-ring $Grid")
    println(io, "│├ snow_depth:               $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       $nlat-ring $Grid")
    println(io, "│└ evaporative_flux:         $nlat-ring $Grid")
    println(io, "├ tracers: $(length(tracer_names)), $tracer_names")
    println(io, "├ particles: $nparticles-element $(typeof(progn.particles))")
    println(io, "├ scale: $(progn.scale[])")
    print(io,   "└ clock: $(progn.clock.time)")
end

"""$(TYPEDSIGNATURES)
Copies entries of `progn_old` into `progn_new`."""
function Base.copy!(progn_new::PrognosticVariables{NF,AT,ST,NSTEPS}, progn_old::PrognosticVariables{NF,AT,ST,NSTEPS}) where {NF,AT,ST,NSTEPS}

    # Core variables using broadcast
    @inbounds for i in 1:NSTEPS
        progn_new.vor[i] .= progn_old.vor[i]
        progn_new.div[i] .= progn_old.div[i]
        progn_new.temp[i] .= progn_old.temp[i]
        progn_new.humid[i] .= progn_old.humid[i]
        progn_new.pres[i] .= progn_old.pres[i]
    end

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

function Base.zero(progn::PrognosticVariables{NF, ArrayType, SpectrumType, nsteps, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector}) where {NF, ArrayType, SpectrumType, nsteps, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ParticleVector}
    (; spectrum, nlat_half, nlayers, nlayers_soil, nparticles) = progn
    
    # initialize regular progn variables 
    progn_new = PrognosticVariables{NF, ArrayType, SpectrumType, nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, 
        ParticleVector}(;
            spectrum, nlat_half, nlayers, nlayers_soil, nparticles,
        )

    # add tracers with zero 
    for (key, value) in progn.tracers 
        progn_new.tracers[key] = ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)
    end 

    # use the same scale 
    progn_new.scale[] = progn.scale[]

    return progn_new
end 

function Base.fill!(progn::PrognosticVariables{NF}, value::Number) where NF

    value_NF = NF(value)

    for i in eachindex(progn.vor)   # each leapfrog time step
        progn.vor[i] .= value_NF
        progn.div[i] .= value_NF
        progn.temp[i] .= value_NF
        progn.humid[i] .= value_NF
        progn.pres[i] .= value_NF
    end

    # ocean
    progn.ocean.sea_surface_temperature .= value_NF
    progn.ocean.sea_ice_concentration .= value_NF
    progn.ocean.sensible_heat_flux .= value_NF
    progn.ocean.evaporative_flux .= value_NF

    # land
    progn.land.soil_temperature .= value_NF
    progn.land.snow_depth .= value_NF
    progn.land.soil_moisture .= value_NF
    progn.land.sensible_heat_flux .= value_NF
    progn.land.evaporative_flux .= value_NF

    # fill tracers
    for (key, value) in progn.tracers 
        for value_i in value # istep of nsteps tuple 
            value_i .= value_NF
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
    progn::PrognosticVariables{NF, ArrayType, SpectrumType, nsteps, SpectralVariable2D, SpectralVariable3D},
    tracers::Tracer...
) where {
        NF,                     # number format
        ArrayType,
        SpectrumType,
        nsteps,
        SpectralVariable2D,
        SpectralVariable3D,
    }
    (; spectrum, nlayers) = progn
    for tracer in tracers
        progn.tracers[tracer.name] = ntuple(i -> zeros(SpectralVariable3D, spectrum, nlayers), nsteps)
    end
end

"""$(TYPEDSIGNATURES)
Delete `tracers` in the prognostic variables `progn` in `progn.tracers::Dict`."""
function Base.delete!(progn::PrognosticVariables, tracers::Tracer...)
    for tracer in tracers
        delete!(progn.tracers, tracer.name)
    end
end