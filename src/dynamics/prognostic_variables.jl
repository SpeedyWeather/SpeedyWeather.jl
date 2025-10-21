function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

export PrognosticVariablesOcean
@kwdef struct PrognosticVariablesOcean{
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

    "Prescribed ocean humidity flux [kg/s/m²]"
    surface_humidity_flux::GridVariable2D = zeros(GridVariable2D, grid)
end

Adapt.@adapt_structure PrognosticVariablesOcean

export PrognosticVariablesLand
@kwdef struct PrognosticVariablesLand{
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

    "Prescribed land sensible heat flux [W/m²], positive up, zero if not used"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Prescribed land humidity flux [kg/s/m²], positive up, zero if not used"
    surface_humidity_flux::GridVariable2D = zeros(GridVariable2D, grid)
end

Adapt.@adapt_structure PrognosticVariablesLand
Base.eltype(::PrognosticVariablesLand{GridType, GridVariable}) where {GridType, GridVariable <: AbstractField{NF}} where {NF} = NF

export PrognosticVariables
@kwdef struct PrognosticVariables{
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <: AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    SpectralVariable4D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    TracerTuple,            # <: NamedTuple{Tuple{Symbol}, Tuple{SpectralVariable4D}}
    ParticleVector,         # <: AbstractVector{Particle{NF}}
    RefValueNF,               # <: Base.RefValue{NF}
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
    ocean::PrognosticVariablesOcean{GridType, GridVariable2D} =
        PrognosticVariablesOcean{GridType, GridVariable2D}(; grid)
    
    "Land variables, soil temperature, snow, and soil moisture"
    land::PrognosticVariablesLand{GridType, GridVariable2D, GridVariable3D} =
        PrognosticVariablesLand{GridType, GridVariable2D, GridVariable3D}(; grid, nlayers=nlayers_soil)

    "Tracers, last dimension is for n tracers [?]"
    tracers::TracerTuple = NamedTuple()

    "Particles for particle advection"
    particles::ParticleVector = zeros(ParticleVector, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::RefValueNF = Ref(one(eltype(land)))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::Clock = Clock()
end

Adapt.@adapt_structure PrognosticVariables

Base.eltype(::PrognosticVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D}) where {SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D <: AbstractField{NF}} where NF = NF 
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
function PrognosticVariables(SG::SpectralGrid{Architecture, SpectrumType, GridType}; tracers::TRACER_DICT=TRACER_DICT(), nsteps=DEFAULT_NSTEPS) where {Architecture, SpectrumType, GridType}

    (; NF, spectrum, grid, nlayers, nlayers_soil, nparticles) = SG
    (; SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector) = SG
    
    tracer_tuple = (; [key => zeros(SpectralVariable4D, spectrum, nlayers, nsteps) for key in keys(tracers)]...)

    return PrognosticVariables{SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, typeof(tracer_tuple), ParticleVector, Base.RefValue{NF}}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps, tracer_tuple
        )
end

"""$(TYPEDSIGNATURES)
Generator function."""
PrognosticVariables(SG::SpectralGrid, model::AbstractModel) = PrognosticVariables(SG, tracers = model.tracers, nsteps = model.time_stepping.nsteps)

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
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sea_ice_concentration:    Field{$NF} on $nlat-ring $Grid")
    println(io, "│├ sensible_heat_flux:       Field{$NF} on $nlat-ring $Grid")
    println(io, "│└ surface_humidity_flux:    Field{$NF} on $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
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
    progn_new.ocean.surface_humidity_flux .= progn_old.ocean.surface_humidity_flux

    # Land variables
    progn_new.land.soil_temperature .= progn_old.land.soil_temperature
    progn_new.land.snow_depth .= progn_old.land.snow_depth
    progn_new.land.soil_moisture .= progn_old.land.soil_moisture
    progn_new.land.sensible_heat_flux .= progn_old.land.sensible_heat_flux
    progn_new.land.surface_humidity_flux .= progn_old.land.surface_humidity_flux

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
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector,
        }) where {
        SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, ParticleVector,
        }

    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn
    
    # initialize regular progn variables 
    progn_new = PrognosticVariables{SpectrumType, GridType,
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, 
        ParticleVector, typeof(progn.scale)}(;
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

"""$(TYPEDSIGNATURES)
Add `tracers` to the prognostic variables `progn` in `progn.tracers::NamedTuple`, reconstructs the `PrognosticVariables`."""
function add(
    progn::PrognosticVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, NF},
    tracers::Tracer...
    ) where {SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, NF}

    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn
    
    new_tracers = union(keys(progn.tracers), keys(tracers))

    progn_new = PrognosticVariables{SpectrumType, GridType, 
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, 
        GridVariable2D, GridVariable3D, typeof(new_tracers), ParticleVector, NF}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps, new_tracers
        )

    copy!(progn_new, progn)
    return progn_new 
end

"""$(TYPEDSIGNATURES)
Delete `tracers` in the prognostic variables `progn` in `progn.tracers::NamedTuple`, reconstructs the `PrognosticVariables`."""
function delete(progn::PrognosticVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, NF}, 
    tracers::Tracer...
    ) where {SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, GridVariable2D, GridVariable3D, TracerTuple, ParticleVector, NF}

    (; spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps) = progn
    
    new_tracers = setdiff(keys(progn.tracers), keys(tracers))
    
    progn_new = PrognosticVariables{SpectrumType, GridType, 
        SpectralVariable2D, SpectralVariable3D, SpectralVariable4D, 
        GridVariable2D, GridVariable3D, typeof(new_tracers), ParticleVector, NF}(;
            spectrum, grid, nlayers, nlayers_soil, nparticles, nsteps, new_tracers
        )

    copy!(progn_new, progn)
    return progn_new
end