function Base.show(io::IO, A::AbstractDiagnosticVariables)
    println(io, "$(typeof(A).name.wrapper)")
    keys = propertynames(A)
    keys_filtered = filter(key -> ~(getfield(A, key) isa Integer), keys)
    n = length(keys_filtered)
    for (i, key) in enumerate(keys_filtered)
        last = i == n
        val = getfield(A, key)
        T = typeof(val)

        if T <: AbstractGridArray
            NF = first_parameter(T)
            nlat = RingGrids.get_nlat(val)
            Grid = RingGrids.nonparametric_type(T)
            s = Base.dims2string(size(val))*", $nlat-ring $Grid{$NF}"
        elseif T <: LowerTriangularArray
            NF = first_parameter(T)
            trunc = val.n - 1
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
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
} <: AbstractDiagnosticVariables

    trunc::Int              # spectral resolution: maximum degree and order of spherical harmonics
    nlat_half::Int          # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    # SPECTRAL TENDENCIES
    "Vorticity of horizontal wind field [1/s]"
    vor_tend  ::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Divergence of horizontal wind field [1/s]"
    div_tend  ::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Absolute temperature [K]"
    temp_tend ::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Specific humidity [kg/kg]"
    humid_tend::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Zonal velocity [m/s]"
    u_tend    ::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Meridional velocity [m/s]"
    v_tend    ::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    "Logarithm of surface pressure [Pa]"
    pres_tend ::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)
    "Tracers [?]"
    tracers_tend::Dict{Symbol, SpectralVariable3D} = Dict{Symbol, SpectralVariable3D}()

    "Zonal velocity [m/s], grid"
    u_tend_grid     ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Meridinoal velocity [m/s], grid"
    v_tend_grid     ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Absolute temperature [K], grid"
    temp_tend_grid  ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Specific humidity [kg/kg], grid"
    humid_tend_grid ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Logarith of surface pressure [Pa], grid"
    pres_tend_grid  ::GridVariable2D = zeros(GridVariable2D, nlat_half)
    "Tracers [?], grid"
    tracers_tend_grid::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function Tendencies(SG::SpectralGrid)
    (; trunc, nlat_half, nlayers, NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return Tendencies{
        NF, ArrayType,
        SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D,
    }(;
        trunc, nlat_half, nlayers
    )
end

export GridVariables

"""Transformed prognostic variables (and u, v, temp_virt) into grid-point space.
$TYPEDFIELDS."""
@kwdef struct GridVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
} <: AbstractDiagnosticVariables

    nlat_half::Int          # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    "Relative vorticity of the horizontal wind [1/s]"
    vor_grid        ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Divergence of the horizontal wind [1/s]"
    div_grid        ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Absolute temperature [K]"
    temp_grid       ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Virtual tempereature [K]"
    temp_virt_grid  ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Specific_humidity [kg/kg]"
    humid_grid      ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Zonal velocity [m/s]"
    u_grid          ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Meridional velocity [m/s]"
    v_grid          ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Logarithm of surface pressure [Pa]"
    pres_grid       ::GridVariable2D = zeros(GridVariable2D, nlat_half)
    "Tracers [?]"
    tracers_grid    ::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()

    "Random pattern controlled by random process [1]"
    random_pattern  ::GridVariable2D = zeros(GridVariable3D, nlat_half)

    # PREVIOUS TIME STEP
    "Absolute temperature [K] at previous time step"
    temp_grid_prev  ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Specific humidity [kg/kg] at previous time step"
    humid_grid_prev ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Zonal velocity [m/s] at previous time step"
    u_grid_prev     ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Meridional velocity [m/s] at previous time step"
    v_grid_prev     ::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    "Logarithm of surface pressure [Pa] at previous time step"
    pres_grid_prev  ::GridVariable2D = zeros(GridVariable2D, nlat_half)
    "Tracers [?] at previous time step"
    tracers_grid_prev::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function GridVariables(SG::SpectralGrid)
    (; nlat_half, nlayers, NF, ArrayType) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return GridVariables{NF, ArrayType, GridVariable2D, GridVariable3D}(;
            nlat_half, nlayers,
        )
end

export DynamicsVariables

"""Intermediate quantities for the dynamics of a given layer.
$(TYPEDFIELDS)"""
@kwdef struct DynamicsVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
} <: AbstractDiagnosticVariables
    
    trunc::Int              # spectral resolution: maximum degree and order of spherical harmonics
    nlat_half::Int          # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    "Multi-purpose a, 3D work array to be reused in various places"
    a::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    
    "Multi-purpose b, 3D work array to be reused in various places"
    b::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
    
    "Multi-purpose a, 3D work array to be reused in various places"
    a_grid::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
    "Multi-purpose b, 3D work array to be reused in various places"
    b_grid::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
    "Multi-purpose a, work array to be reused in various places"
    a_2D::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)
    
    "Multi-purpose b, work array to be reused in various places"
    b_2D::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)
    
    "Multi-purpose a, work array to be reused in various places"
    a_2D_grid::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Multi-purpose b, work array to be reused in various places"
    b_2D_grid::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Pressure flux (uₖ, vₖ)⋅∇ln(pₛ)"
    uv∇lnp::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
    "Sum of Δσₖ-weighted uv∇lnp above"
    uv∇lnp_sum_above::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
    "Sum of div_weighted from top to k"
    div_sum_above::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)
    
    "Virtual temperature [K], spectral for geopotential"
    temp_virt::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)

    "Geopotential [m²/s²] on full layers"
    geopot::SpectralVariable3D = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)

    "Vertical velocity (dσ/dt), on half levels k+1/2 below, pointing to the surface (σ=1)"
    σ_tend::GridVariable3D = zeros(GridVariable3D, nlat_half, nlayers)

    "Zonal gradient of log surf pressure"
    ∇lnp_x::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Meridional gradient of log surf pressure"
    ∇lnp_y::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Vertical average of zonal velocity [m/s]"
    u_mean_grid::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Vertical average of meridional velocity [m/s]"
    v_mean_grid::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Vertical average of divergence [1/s], grid"
    div_mean_grid::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Vertical average of divergence [1/s], spectral"
    div_mean::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)

    "Scratch memory for the transforms"
    scratch_memory = SpeedyTransforms.ScratchMemory(NF, ArrayType, nlat_half, GridVariable2D, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator function."""
function DynamicsVariables(SG::SpectralGrid)
    (; trunc, nlat_half, nlayers, NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return DynamicsVariables{NF, ArrayType, SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D}(;
            trunc, nlat_half, nlayers,
        )
end

export PhysicsVariables

"""
Diagnostic variables of the physical parameterizations.
$(TYPEDFIELDS)"""
@kwdef struct PhysicsVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
} <: AbstractDiagnosticVariables

    nlat_half::Int

    # PRECIPITATION
    "Accumulated large-scale precipitation [m]"
    precip_large_scale::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Accumulated large-scale precipitation [m]"
    precip_convection::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Rate of large-scale precipitation [m/s], instantaneous"
    precip_rate_large_scale::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Rate of large-scale precipitation [m/s], instantaneous"
    precip_rate_convection::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Cloud top [m]"
    cloud_top::GridVariable2D = zeros(GridVariable2D, nlat_half)            
    
    # LAND
    "Availability of soil moisture to evaporation [1]"
    soil_moisture_availability::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "River runoff [m/s], diagnostic overflow from soil moisture"
    river_runoff::GridVariable2D = zeros(GridVariable2D, nlat_half)

    # SURFACE FLUXES
    "Sensible heat flux [W/m²], positive up"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)
    sensible_heat_flux_ocean::GridVariable2D = zeros(GridVariable2D, nlat_half)
    sensible_heat_flux_land::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Evaporative flux [kg/s/m^2], positive up"
    evaporative_flux::GridVariable2D = zeros(GridVariable2D, nlat_half)
    evaporative_flux_ocean::GridVariable2D = zeros(GridVariable2D, nlat_half)
    evaporative_flux_land::GridVariable2D = zeros(GridVariable2D, nlat_half)

    # RADIATION
    "Surface radiation: shortwave up [W/m²]"
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Surface radiation: shortwave down [W/m²]"
    surface_shortwave_down::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Surface radiation: longwave up [W/m²]"
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, nlat_half)
    surface_longwave_up_ocean::GridVariable2D = zeros(GridVariable2D, nlat_half)
    surface_longwave_up_land::GridVariable2D = zeros(GridVariable2D, nlat_half)
    
    "Surface radiation: longwave down [W/m²]"
    surface_longwave_down::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Outgoing shortwave radiation [W/m^2]"
    outgoing_shortwave_radiation::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Outgoing longwave radiation [W/m^2]"
    outgoing_longwave_radiation::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Cosine of solar zenith angle [1]"
    cos_zenith::GridVariable2D = zeros(GridVariable2D, nlat_half)           
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PhysicsVariables(SG::SpectralGrid)
    (; nlat_half, NF, ArrayType) = SG
    (; GridVariable2D) = SG

    return PhysicsVariables{NF, ArrayType, GridVariable2D}(; nlat_half)
end

export ParticleVariables

"""Diagnostic variables for the particle advection
$(TYPEDFIELDS)"""
@kwdef struct ParticleVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    ParticleVector,         # <: AbstractGridArray
    VectorNF,               # Vector{NF} or CuVector{NF}
    Grid,                   # <:AbstractGridArray
} <: AbstractDiagnosticVariables
    "Number of particles"
    nparticles::Int

    "Number of latitudes on one hemisphere (Eq. incld.), resolution parameter of Grid"
    nlat_half::Int

    "Work array: particle locations"
    locations::ParticleVector = zeros(ParticleVector, nparticles)

    "Work array: velocity u"
    u::VectorNF = VectorNF(zeros(nparticles))

    "Work array: velocity v"
    v::VectorNF = VectorNF(zeros(nparticles))

    "Work array: velocity w = dσ/dt"
    σ_tend::VectorNF = VectorNF(zeros(nparticles))

    "Interpolator to interpolate velocity fields onto particle positions"
    interpolator::AnvilInterpolator{NF, Grid} = AnvilInterpolator(NF, Grid, nlat_half, nparticles)
end

"""$(TYPEDSIGNATURES)
Generator function."""
function ParticleVariables(SG::SpectralGrid)
    (; nlat_half, nparticles, NF, ArrayType, Grid) = SG
    (; ParticleVector) = SG
    VectorNF = ArrayType{NF, 1}
    return ParticleVariables{NF, ArrayType, ParticleVector, VectorNF, Grid}(; nlat_half, nparticles)
end

export DiagnosticVariables

"""All diagnostic variables.
$(TYPEDFIELDS)"""
struct DiagnosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    Grid,                   # <:AbstractGridArray
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
    ParticleVector,         # <: AbstractGridArray
    VectorType,             # <: AbstractVector
    MatrixType,             # <: AbstractMatrix
} <: AbstractDiagnosticVariables

    # DIMENSIONS
    "Spectral resolution: Max degree of spherical harmonics (0-based)"
    trunc::Int

    "Grid resoltion: Number of latitude rings on one hemisphere (Equator incl.)"
    nlat_half::Int

    "Number of vertical layers"
    nlayers::Int

    "Number of particles for particle advection"
    nparticles::Int

    "Tendencies (spectral and grid) of the prognostic variables"
    tendencies::Tendencies{NF, ArrayType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D}
    
    "Gridded prognostic variables"
    grid::GridVariables{NF, ArrayType, GridVariable2D, GridVariable3D}
    
    "Intermediate variables for the dynamical core"
    dynamics::DynamicsVariables{NF, ArrayType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D}
    
    "Global fields returned from physics parameterizations"
    physics::PhysicsVariables{NF, ArrayType, GridVariable2D}
    
    "Intermediate variables for the particle advection"
    particles::ParticleVariables{NF, ArrayType, ParticleVector, VectorType, Grid}
    
    "Vertical column for the physics parameterizations"
    column::ColumnVariables{NF, VectorType, MatrixType}

    "Average temperature of every horizontal layer [K]"
    temp_average::VectorType

    "Scale applied to vorticity and divergence"
    scale::Base.RefValue{NF}
end

function DiagnosticVariables(SG::SpectralGrid, model::Union{Barotropic, ShallowWater})
    diagn = DiagnosticVariables(SG)
    add!(diagn, model.tracers)
    return diagn
end

# decide on spectral resolution `nbands` of radiation schemes
function DiagnosticVariables(SG::SpectralGrid, model::PrimitiveEquation)
    nbands_shortwave = get_nbands(model.shortwave_radiation)
    nbands_longwave = get_nbands(model.longwave_radiation)
    diagn =  DiagnosticVariables(SG; nbands_shortwave, nbands_longwave)
    add!(diagn, model.tracers)
    return diagn
end

"""$(TYPEDSIGNATURES)
Generator function."""
function DiagnosticVariables(
    SG::SpectralGrid;
    nbands_shortwave::Integer = 0,
    nbands_longwave::Integer = 0,
)
    (; trunc, nlat_half, nparticles, NF, nlayers) = SG

    tendencies = Tendencies(SG)
    grid = GridVariables(SG)
    dynamics = DynamicsVariables(SG)
    physics = PhysicsVariables(SG)
    particles = ParticleVariables(SG)
    column = ColumnVariables(SG; nbands_shortwave, nbands_longwave)
    temp_average = SG.VectorType(undef, nlayers)

    scale = Ref(one(NF))

    return DiagnosticVariables(
        trunc, nlat_half, nlayers, nparticles,
        tendencies, grid, dynamics, physics, particles,
        column, temp_average, scale,
    )
end

function Base.show(
    io::IO,
    diagn::DiagnosticVariables{NF, ArrayType, Grid},
) where {NF, ArrayType, Grid}
    println(io, "DiagnosticVariables{$NF, $ArrayType, $Grid}")
    
    (; trunc, nlat_half, nlayers, nparticles) = diagn
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    ntracers = length(diagn.grid.tracers_grid)
    println(io, "├ resolution: T$trunc, $nlat rings, $nlayers layers, $ntracers tracers, $nparticles particles")
    println(io, "├ tendencies::Tendencies")
    println(io, "├ grid::GridVariables")
    println(io, "├ dynamics::DynamicsVariables")
    println(io, "├ physics::PhysicsVariables")
    println(io, "├ particles::ParticleVariables")
    println(io, "├ columns::ColumnVariables")
    println(io, "├ temp_average::$(typeof(diagn.temp_average))")
    print(io,   "└ scale: $(diagn.scale[])")
end

function add!(diagn::DiagnosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    Grid,                   # <:AbstractGridArray
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    GridVariable3D,         # <: AbstractGridArray
},
    tracers::Tracer...,
    ) where {
        NF,                     # <: AbstractFloat
        ArrayType,              # Array, CuArray, ...
        Grid,                   # <:AbstractGridArray
        SpectralVariable2D,     # <: LowerTriangularArray
        SpectralVariable3D,     # <: LowerTriangularArray
        GridVariable2D,         # <: AbstractGridArray
        GridVariable3D,         # <: AbstractGridArray
    }
    (; trunc, nlat_half, nlayers) = diagn
    for tracer in tracers
        diagn.tendencies.tracers_tend[tracer.name] = zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers)
        diagn.tendencies.tracers_tend_grid[tracer.name] = zeros(GridVariable3D, nlat_half, nlayers)
        diagn.grid.tracers_grid[tracer.name] = zeros(GridVariable3D, nlat_half, nlayers)
        diagn.grid.tracers_grid_prev[tracer.name] = zeros(GridVariable3D, nlat_half, nlayers)
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

RingGrids.eachgridpoint(diagn::DiagnosticVariables) = eachgridpoint(diagn.grid.vor_grid)