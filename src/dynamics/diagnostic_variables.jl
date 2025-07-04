function Base.show(io::IO, A::AbstractDiagnosticVariables)
    println(io, "$(typeof(A).name.wrapper)")
    keys = propertynames(A)
    keys_filtered = filter(key -> ~(getfield(A, key) isa Integer), keys)
    n = length(keys_filtered)
    for (i, key) in enumerate(keys_filtered)
        last = i == n
        val = getfield(A, key)
        T = typeof(val)

        if T <: AbstractField
            NF = first_parameter(T)
            nlat = RingGrids.get_nlat(val)
            Grid = RingGrids.nonparametric_type(T)
            s = Base.dims2string(size(val))*", $nlat-ring $Grid{$NF}"
        elseif T <: LowerTriangularArray
            NF = first_parameter(T)
            trunc = LowerTriangularArrays.truncation(val.spectrum)
            s = Base.dims2string(size(val, as=Matrix))*", T$trunc LowerTriangularArray{$NF}"
        else
            s = "$T"
        end

        ~last ? println(io, "├ $key: $s") :
                print(io,  "└ $key: $s")
    end
end

first_parameter(::Type{<:AbstractArray{T}}) where T = T

export Tendencies

"""Tendencies of the prognostic variables in spectral and grid-point space
$(TYPEDFIELDS)"""
@kwdef struct Tendencies{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <: AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
} <: AbstractDiagnosticVariables

    spectrum::SpectrumType            # spectral resolution: maximum degree and order of spherical harmonics
    grid::GridType              # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    # SPECTRAL TENDENCIES
    "Vorticity of horizontal wind field [1/s]"
    vor_tend  ::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Divergence of horizontal wind field [1/s]"
    div_tend  ::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Absolute temperature [K]"
    temp_tend ::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Specific humidity [kg/kg]"
    humid_tend::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Zonal velocity [m/s]"
    u_tend    ::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Meridional velocity [m/s]"
    v_tend    ::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    "Logarithm of surface pressure [Pa]"
    pres_tend ::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)
    "Tracers [?]"
    tracers_tend::Dict{Symbol, SpectralVariable3D} = Dict{Symbol, SpectralVariable3D}()

    "Zonal velocity [m/s], grid"
    u_tend_grid     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Meridinoal velocity [m/s], grid"
    v_tend_grid     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Absolute temperature [K], grid"
    temp_tend_grid  ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Specific humidity [kg/kg], grid"
    humid_tend_grid ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Logarith of surface pressure [Pa], grid"
    pres_tend_grid  ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?], grid"
    tracers_tend_grid::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function Tendencies(SG::SpectralGrid)
    (; spectrum, grid, nlayers, NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return Tendencies{
        NF, ArrayType, typeof(spectrum), typeof(grid),
        SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D,
    }(;
        spectrum, grid, nlayers
    )
end

export GridVariables

"""Transformed prognostic variables (and u, v, temp_virt) into grid-point space.
$TYPEDFIELDS."""
@kwdef struct GridVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridType,               # <:AbstractGrid
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
} <: AbstractDiagnosticVariables

    grid::GridType              # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    "Relative vorticity of the horizontal wind [1/s]"
    vor_grid        ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Divergence of the horizontal wind [1/s]"
    div_grid        ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Absolute temperature [K]"
    temp_grid       ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Virtual tempereature [K]"
    temp_virt_grid  ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Specific_humidity [kg/kg]"
    humid_grid      ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Zonal velocity [m/s]"
    u_grid          ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Meridional velocity [m/s]"
    v_grid          ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Logarithm of surface pressure [Pa]"
    pres_grid       ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?]"
    tracers_grid    ::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()

    "Random pattern controlled by random process [1]"
    random_pattern  ::GridVariable2D = zeros(GridVariable2D, grid)

    # PREVIOUS TIME STEP
    "Absolute temperature [K] at previous time step"
    temp_grid_prev  ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Specific humidity [kg/kg] at previous time step"
    humid_grid_prev ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Zonal velocity [m/s] at previous time step"
    u_grid_prev     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Meridional velocity [m/s] at previous time step"
    v_grid_prev     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Logarithm of surface pressure [Pa] at previous time step"
    pres_grid_prev  ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?] at previous time step"
    tracers_grid_prev::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function GridVariables(SG::SpectralGrid)
    (; grid, nlayers, NF, ArrayType) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return GridVariables{NF, ArrayType, typeof(grid), GridVariable2D, GridVariable3D}(;
            grid, nlayers,
        )
end

export DynamicsVariables

"""Intermediate quantities for the dynamics of a given layer.
$(TYPEDFIELDS)"""
@kwdef struct DynamicsVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    ScratchMemoryType,      # <: ScratchMemory{NF, ArrayType{Complex{NF},3}}
} <: AbstractDiagnosticVariables
    
    spectrum::SpectrumType            # spectral resolution: maximum degree and order of spherical harmonics
    grid::GridType              # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    "Multi-purpose a, 3D work array to be reused in various places"
    a::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    
    "Multi-purpose b, 3D work array to be reused in various places"
    b::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)
    
    "Multi-purpose a, 3D work array to be reused in various places"
    a_grid::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Multi-purpose b, 3D work array to be reused in various places"
    b_grid::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Multi-purpose a, work array to be reused in various places"
    a_2D::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)
    
    "Multi-purpose b, work array to be reused in various places"
    b_2D::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)
    
    "Multi-purpose a, work array to be reused in various places"
    a_2D_grid::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Multi-purpose b, work array to be reused in various places"
    b_2D_grid::GridVariable2D = zeros(GridVariable2D, grid)

    "Pressure flux (uₖ, vₖ)⋅∇ln(pₛ)"
    uv∇lnp::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Sum of Δσₖ-weighted uv∇lnp above"
    uv∇lnp_sum_above::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Sum of div_weighted from top to k"
    div_sum_above::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    
    "Virtual temperature [K], spectral for geopotential"
    temp_virt::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)

    "Geopotential [m²/s²] on full layers"
    geopot::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)

    "Vertical velocity (dσ/dt), on half levels k+1/2 below, pointing to the surface (σ=1)"
    σ_tend::GridVariable3D = zeros(GridVariable3D, grid, nlayers)

    "Zonal gradient of log surf pressure"
    ∇lnp_x::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Meridional gradient of log surf pressure"
    ∇lnp_y::GridVariable2D = zeros(GridVariable2D, grid)

    "Vertical average of zonal velocity [m/s]"
    u_mean_grid::GridVariable2D = zeros(GridVariable2D, grid)

    "Vertical average of meridional velocity [m/s]"
    v_mean_grid::GridVariable2D = zeros(GridVariable2D, grid)

    "Vertical average of divergence [1/s], grid"
    div_mean_grid::GridVariable2D = zeros(GridVariable2D, grid)

    "Vertical average of divergence [1/s], spectral"
    div_mean::SpectralVariable2D = zeros(SpectralVariable2D, spectrum)

    "Scratch memory for the transforms"
    scratch_memory::ScratchMemoryType = SpeedyTransforms.ScratchMemory(NF, ArrayType, grid, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator function. If a `spectral_transform` is handed over, the same scratch memory is used."""
function DynamicsVariables(SG::SpectralGrid; 
                           spectral_transform::Union{Nothing,SpectralTransform}=nothing) 
    (; spectrum, grid, nlayers, NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    if isnothing(spectral_transform)
        return DynamicsVariables{NF, ArrayType, typeof(spectrum), typeof(grid),
            SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, SpeedyTransforms.ScratchMemory{NF, ArrayType{Complex{NF}, 3}}}(;
                spectrum, grid, nlayers,
            )
    else 
        scratch_memory = spectral_transform.scratch_memory 

        return DynamicsVariables{NF, ArrayType, typeof(spectrum), typeof(grid),
            SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, typeof(scratch_memory)}(;
                spectrum, grid, nlayers, scratch_memory
            )
    end 
end


export DynamicsVariablesOcean
@kwdef struct DynamicsVariablesOcean{
    NF,
    ArrayType,
    GridType,
    GridVariable2D,
} <: AbstractDiagnosticVariables

    grid::GridType
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, grid)
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    albedo::GridVariable2D = zeros(GridVariable2D, grid)
end

DynamicsVariablesOcean(SG::SpectralGrid) =
    DynamicsVariablesOcean{SG.NF, SG.ArrayType, typeof(SG.grid), SG.GridVariable2D}(; SG.grid)

export DynamicsVariablesLand
@kwdef struct DynamicsVariablesLand{
    NF,
    ArrayType,
    GridType,
    GridVariable2D,
} <: AbstractDiagnosticVariables

    grid::GridType
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, grid)
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    albedo::GridVariable2D = zeros(GridVariable2D, grid)

    "Availability of soil moisture to evaporation [1]"
    soil_moisture_availability::GridVariable2D = zeros(GridVariable2D, grid)

    "River runoff [m/s], diagnostic overflow from soil moisture"
    river_runoff::GridVariable2D = zeros(GridVariable2D, grid)
end

DynamicsVariablesLand(SG::SpectralGrid) =
    DynamicsVariablesLand{SG.NF, SG.ArrayType, typeof(SG.grid), SG.GridVariable2D}(; SG.grid)


export PhysicsVariables

"""
Diagnostic variables of the physical parameterizations.
$(TYPEDFIELDS)"""
@kwdef struct PhysicsVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridType,               # <:AbstractGrid
    GridVariable2D,         # <: AbstractField
} <: AbstractDiagnosticVariables

    grid::GridType          # resolution of grid

    ocean::DynamicsVariablesOcean{NF, ArrayType, GridType, GridVariable2D}
    land::DynamicsVariablesLand{NF, ArrayType, GridType, GridVariable2D}

    # PRECIPITATION
    "Accumulated large-scale precipitation [m]"
    precip_large_scale::GridVariable2D = zeros(GridVariable2D, grid)

    "Accumulated large-scale precipitation [m]"
    precip_convection::GridVariable2D = zeros(GridVariable2D, grid)

    "Rate of large-scale precipitation [m/s], instantaneous"
    precip_rate_large_scale::GridVariable2D = zeros(GridVariable2D, grid)

    "Rate of large-scale precipitation [m/s], instantaneous"
    precip_rate_convection::GridVariable2D = zeros(GridVariable2D, grid)

    "Cloud top [m]"
    cloud_top::GridVariable2D = zeros(GridVariable2D, grid)            
    
    # SURFACE FLUXES
    "Sensible heat flux [W/m²], positive up"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Evaporative flux [kg/s/m^2], positive up"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, grid)

    # RADIATION
    "Surface radiation: shortwave up [W/m²]"
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface radiation: shortwave down [W/m²]"
    surface_shortwave_down::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface radiation: longwave up [W/m²]"
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface radiation: longwave down [W/m²]"
    surface_longwave_down::GridVariable2D = zeros(GridVariable2D, grid)

    "Outgoing shortwave radiation [W/m^2]"
    outgoing_shortwave_radiation::GridVariable2D = zeros(GridVariable2D, grid)

    "Outgoing longwave radiation [W/m^2]"
    outgoing_longwave_radiation::GridVariable2D = zeros(GridVariable2D, grid)

    "Albedo [1]"
    albedo::GridVariable2D = zeros(GridVariable2D, grid)

    "Cosine of solar zenith angle [1]"
    cos_zenith::GridVariable2D = zeros(GridVariable2D, grid)           
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PhysicsVariables(SG::SpectralGrid)
    (; grid, NF, ArrayType) = SG
    (; GridVariable2D) = SG

    ocean = DynamicsVariablesOcean(SG)
    land = DynamicsVariablesLand(SG)

    return PhysicsVariables{NF, ArrayType, typeof(grid), GridVariable2D}(; grid, ocean, land)
end

export ParticleVariables

"""Diagnostic variables for the particle advection
$(TYPEDFIELDS)"""
@kwdef struct ParticleVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    ParticleVector,         # <: AbstractField
    VectorNF,               # Vector{NF} or CuVector{NF}
    Interpolator,           # <:AbstractInterpolator
} <: AbstractDiagnosticVariables
    "Number of particles"
    nparticles::Int

    "Work array: particle locations"
    locations::ParticleVector = zeros(ParticleVector, nparticles)

    "Work array: velocity u"
    u::VectorNF = VectorNF(zeros(nparticles))

    "Work array: velocity v"
    v::VectorNF = VectorNF(zeros(nparticles))

    "Work array: velocity w = dσ/dt"
    σ_tend::VectorNF = VectorNF(zeros(nparticles))

    "Interpolator to interpolate velocity fields onto particle positions"
    interpolator::Interpolator
end

"""$(TYPEDSIGNATURES)
Generator function."""
function ParticleVariables(SG::SpectralGrid)
    (; nparticles, NF, ArrayType) = SG
    (; ParticleVector) = SG
    VectorNF = ArrayType{NF, 1}
    interpolator = RingGrids.AnvilInterpolator(SG.grid, nparticles; NF, ArrayType)
    return ParticleVariables{NF,ArrayType, ParticleVector, VectorNF, typeof(interpolator)}(;
            nparticles, interpolator)
end

export DiagnosticVariables

"""All diagnostic variables.
$(TYPEDFIELDS)"""
struct DiagnosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    ParticleVector,         # <: AbstractField
    VectorType,             # <: AbstractVector
    MatrixType,             # <: AbstractMatrix
    ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
    Interpolator,           # <:AbstractInterpolator
} <: AbstractDiagnosticVariables

    # DIMENSIONS
    "Spectral resolution: Max degree of spherical harmonics"
    spectrum::SpectrumType

    "Grid resolution: Number of latitude rings on one hemisphere (Equator incl.)"
    grid_used::GridType          # TODO cannot be named `grid` because of `GridVariables` ...

    "Number of vertical layers"
    nlayers::Int

    "Number of particles for particle advection"
    nparticles::Int

    "Tendencies (spectral and grid) of the prognostic variables"
    tendencies::Tendencies{NF, ArrayType, SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D}
    
    "Gridded prognostic variables"
    grid::GridVariables{NF, ArrayType, GridType, GridVariable2D, GridVariable3D}
    
    "Intermediate variables for the dynamical core"
    dynamics::DynamicsVariables{NF, ArrayType, SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ScratchMemoryType}
    
    "Global fields returned from physics parameterizations"
    physics::PhysicsVariables{NF, ArrayType, GridType, GridVariable2D}
    
    "Intermediate variables for the particle advection"
    particles::ParticleVariables{NF, ArrayType, ParticleVector, VectorType, Interpolator}
    
    "Vertical column for the physics parameterizations"
    column::ColumnVariables{NF, VectorType, MatrixType}

    "Average temperature of every horizontal layer [K]"
    temp_average::VectorType

    "Scale applied to vorticity and divergence"
    scale::Base.RefValue{NF}
end

function DiagnosticVariables(SG::SpectralGrid, model::Union{Barotropic, ShallowWater})
    diagn = DiagnosticVariables(SG; spectral_transform=model.spectral_transform)
    add!(diagn, model.tracers)
    return diagn
end

# decide on spectral resolution `nbands` of radiation schemes
function DiagnosticVariables(SG::SpectralGrid, model::PrimitiveEquation)
    nbands_shortwave = get_nbands(model.shortwave_radiation)
    nbands_longwave = get_nbands(model.longwave_radiation)
    diagn =  DiagnosticVariables(SG; spectral_transform=model.spectral_transform, nbands_shortwave, nbands_longwave)
    add!(diagn, model.tracers)
    return diagn
end

"""$(TYPEDSIGNATURES)
Generator function. If a `transform` is handed over, the same scratch memory is used."""
function DiagnosticVariables(
    SG::SpectralGrid;
    spectral_transform::Union{Nothing, SpectralTransform} = nothing,
    nbands_shortwave::Integer = 0,
    nbands_longwave::Integer = 0,
)
    (; spectrum, nparticles, NF, nlayers) = SG

    tendencies = Tendencies(SG)
    grid = GridVariables(SG)
    dynamics = DynamicsVariables(SG; spectral_transform)
    physics = PhysicsVariables(SG)
    particles = ParticleVariables(SG)
    column = ColumnVariables(SG; nbands_shortwave, nbands_longwave)
    temp_average = SG.VectorType(undef, nlayers)

    scale = Ref(one(NF))

    return DiagnosticVariables(
        spectrum, SG.grid, nlayers, nparticles,
        tendencies, grid, dynamics, physics, particles,
        column, temp_average, scale,
    )
end

function Base.show(
    io::IO,
    diagn::DiagnosticVariables{NF, ArrayType},
) where {NF, ArrayType}
    println(io, "DiagnosticVariables{$NF, $ArrayType}")
    
    (; spectrum, nlayers, nparticles) = diagn
    grid = diagn.grid_used   # TODO grid is used by 'GridVariables'
    nlat = RingGrids.get_nlat(grid)
    Grid = RingGrids.nonparametric_type(grid)

    ntracers = length(diagn.grid.tracers_grid)
    println(io, "├ spectrum: T$(truncation(spectrum)), $nlayers layers, $ntracers tracers")
    println(io, "├ grid: $nlat-ring, $nlayers-layer $Grid")
    println(io, "├ tendencies::Tendencies")
    println(io, "├ grid::GridVariables")
    println(io, "├ dynamics::DynamicsVariables")
    println(io, "├ physics::PhysicsVariables")
    println(io, "├ particles::ParticleVariables, $nparticles particles")
    println(io, "├ columns::ColumnVariables")
    println(io, "├ temp_average::$(typeof(diagn.temp_average))")
    print(io,   "└ scale: $(diagn.scale[])")
end

function add!(diagn::DiagnosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
},
    tracers::Tracer...,
    ) where {
        NF,                     # <: AbstractFloat
        ArrayType,              # Array, CuArray, ...
        SpectrumType,           # <: AbstractSpectrum
        GridType,               # <:AbstractGrid
        SpectralVariable2D,     # <: LowerTriangularArray
        SpectralVariable3D,     # <: LowerTriangularArray
        GridVariable2D,         # <: AbstractField
        GridVariable3D,         # <: AbstractField
    }
    (; spectrum, nlayers) = diagn
    grid = diagn.grid_used   # TODO grid is used for 'GridVariables'
    for tracer in tracers
        diagn.tendencies.tracers_tend[tracer.name] = zeros(SpectralVariable3D, spectrum, nlayers)
        diagn.tendencies.tracers_tend_grid[tracer.name] = zeros(GridVariable3D, grid, nlayers)
        diagn.grid.tracers_grid[tracer.name] = zeros(GridVariable3D, grid, nlayers)
        diagn.grid.tracers_grid_prev[tracer.name] = zeros(GridVariable3D, grid, nlayers)
    end
end

function Base.delete!(diagn::DiagnosticVariables, tracers::Tracer...)
    for tracer in tracers
        delete!(diagn.tendencies.tracers_tend, tracer.name)
        delete!(diagn.tendencies.tracers_tend_grid, tracer.name)
        delete!(diagn.grid.tracers_grid, tracer.name)
        delete!(diagn.grid.tracers_grid_prev, tracer.name)
    end
end

"""$(TYPEDSIGNATURES)
Set the tendencies for the barotropic model to `x`."""
function Base.fill!(tendencies::Tendencies, x, ::Type{<:Barotropic})
    fill!(tendencies.u_tend_grid, x)
    fill!(tendencies.v_tend_grid, x)
    fill!(tendencies.vor_tend, x)

    for tracer in values(tendencies.tracers_tend)
        fill!(tracer, x)
    end

    for tracer in values(tendencies.tracers_tend_grid)
        fill!(tracer, x)
    end

    return tendencies
end

"""$(TYPEDSIGNATURES)
Set the tendencies for the shallow-water model to `x`."""
function Base.fill!(tendencies::Tendencies, x, ::Type{<:ShallowWater})
    fill!(tendencies, x, Barotropic)    # all tendencies also in Barotropic
    fill!(tendencies.div_tend, x)       # plus divergence and pressure
    fill!(tendencies.pres_tend_grid, x)
    fill!(tendencies.pres_tend, x)
    return tendencies
end

"""$(TYPEDSIGNATURES)
Set the tendencies for the primitive dry model to `x`."""
function Base.fill!(tendencies::Tendencies, x, ::Type{<:PrimitiveDry})
    fill!(tendencies, x, ShallowWater)  # all tendencies also in ShallowWater
    fill!(tendencies.temp_tend, x)      # plus temperature 
    fill!(tendencies.temp_tend_grid, x)
    return tendencies
end

"""
$(TYPEDSIGNATURES)
Set the tendencies for the primitive wet model to `x`."""
function Base.fill!(tendencies::Tendencies, x, ::Type{<:PrimitiveWet})
    fill!(tendencies, x, PrimitiveDry)  # all tendencies also in PrimitiveDry
    fill!(tendencies.humid_tend, x)     # plus humidity
    fill!(tendencies.humid_tend_grid, x)
    return tendencies
end

# fallback to primitive wet
Base.fill!(tendencies::Tendencies, x) = Base.fill!(tendencies, x, PrimitiveWet)
Base.fill!(tendencies::Tendencies, x, model::AbstractModel) = Base.fill!(tendencies, x, typeof(model))

RingGrids.eachlayer(diagn::DiagnosticVariables) = eachlayer(diagn.grid.vor_grid)
RingGrids.eachgridpoint(diagn::DiagnosticVariables) = eachgridpoint(diagn.grid.vor_grid)
RingGrids.eachindex(diagn::DiagnosticVariables) = eachindex(diagn.grid.vor_grid)