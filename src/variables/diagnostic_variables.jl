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
            Grid = nonparametric_type(T)
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
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <: AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    SpectralTracerTuple,    # <: NamedTuple{Tuple{Symbol}, Tuple{SpectralVariable3D}}
    GridTracerTuple,        # <: NamedTuple{Tuple{Symbol}, Tuple{GridVariable3D}}
} <: AbstractDiagnosticVariables

    spectrum::SpectrumType            # spectral resolution: maximum degree and order of spherical harmonics
    grid::GridType                    # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int                      # number of vertical layers

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
    tracers_tend::SpectralTracerTuple = NamedTuple()

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
    tracers_tend_grid::GridTracerTuple = NamedTuple()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function Tendencies(SG::SpectralGrid; tracers::TRACER_DICT=TRACER_DICT())
    (; spectrum, grid, nlayers) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    tracers_tend = (; [key => zeros(SpectralVariable3D, spectrum, nlayers) for key in keys(tracers)]...)

    tracers_tend_grid = (; [key => zeros(GridVariable3D, grid, nlayers) for key in keys(tracers)]...)

    return Tendencies{
        typeof(spectrum), typeof(grid),
        SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D,
        typeof(tracers_tend), typeof(tracers_tend_grid),
    }(;
        spectrum, grid, nlayers,
        tracers_tend, tracers_tend_grid
    )
end

Adapt.@adapt_structure Tendencies

export GridVariables

"""Transformed prognostic variables (and u, v, temp_virt) into grid-point space.
$TYPEDFIELDS."""
@kwdef struct GridVariables{
    GridType,               # <:AbstractGrid
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    GridTracerTuple,        # <: NamedTuple{Tuple{Symbol}, Tuple{GridVariable3D}}
} <: AbstractDiagnosticVariables

    grid::GridType          # grid resolution: number of latitude rings on one hemisphere (Eq. incl.)
    nlayers::Int            # number of vertical layers

    "Relative vorticity of the horizontal wind [1/s]"
    vor_grid        ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Divergence of the horizontal wind [1/s]"
    div_grid        ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Absolute temperature [K]"
    temp_grid       ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Virtual temperature [K]"
    temp_virt_grid  ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Specific_humidity [kg/kg]"
    humid_grid      ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Zonal velocity [m/s]"
    u_grid          ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Meridional velocity [m/s]"
    v_grid          ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Logarithm of surface pressure [Pa]"
    pres_grid       ::GridVariable2D = zeros(GridVariable2D, grid)
    "Geopotential [m²/s²]"
    geopotential    ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Random pattern controlled by random process [1]"
    random_pattern  ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?]"
    tracers_grid    ::GridTracerTuple = NamedTuple()

    # PREVIOUS TIME STEP
    "Absolute temperature [K] at previous time step"
    temp_grid_prev  ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Specific humidity [kg/kg] at previous time step"
    humid_grid_prev ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Zonal velocity [m/s] at previous time step"
    u_grid_prev     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Meridional velocity [m/s] at previous time step"
    v_grid_prev     ::GridVariable3D = zeros(GridVariable3D, grid, nlayers)
    "Surface pressure [Pa] at previous time step (not logarithm!)"
    pres_grid_prev  ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?] at previous time step"
    tracers_grid_prev::GridTracerTuple = NamedTuple()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function GridVariables(SG::SpectralGrid; tracers::TRACER_DICT=TRACER_DICT())
    (; grid, nlayers) = SG
    (; GridVariable2D, GridVariable3D) = SG

    tracers_grid = (; [key => zeros(GridVariable3D, grid, nlayers) for key in keys(tracers)]...)

    tracers_grid_prev = (; [key => zeros(GridVariable3D, grid, nlayers) for key in keys(tracers)]...)

    return GridVariables{typeof(grid), GridVariable2D, GridVariable3D, typeof(tracers_grid)}(;
            grid, nlayers,
            tracers_grid, tracers_grid_prev
        )
end

Adapt.@adapt_structure GridVariables

export DynamicsVariables

"""Intermediate quantities for the dynamics of a given layer.
$(TYPEDFIELDS)"""
@kwdef struct DynamicsVariables{
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    ScratchMemoryType,      # <: ScratchMemory{ArrayType{Complex{NF},3}}
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
    
    "Virtual temperature [K], spectral, used for geopotential"
    temp_virt::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)

    "Geopotential [m²/s²] on full layers"
    geopotential::SpectralVariable3D = zeros(SpectralVariable3D, spectrum, nlayers)

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
    (; architecture, spectrum, grid, nlayers, NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    if isnothing(spectral_transform)    # then create ScratchMemory now
        scratch_memory = SpeedyTransforms.ScratchMemory(NF, architecture, grid, nlayers)
    else                                # otherwise reuse existing ScratchMemory
        scratch_memory = spectral_transform.scratch_memory
    end

    return DynamicsVariables{typeof(spectrum), typeof(grid),
        SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D, typeof(scratch_memory)}(;
            spectrum, grid, nlayers, scratch_memory
        )
end

Adapt.@adapt_structure DynamicsVariables

export ParticleVariables

"""Diagnostic variables for the particle advection
$(TYPEDFIELDS)"""
@kwdef struct ParticleVariables{
    ParticleVector,         # <: AbstractField
    VectorNF,               # Vector{NF} or CuVector{NF}
    LocatorType,           # <:AbstractLocator
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

    "Locator of the interpolation holding work arrays to interpolate velocity fields onto particle positions"
    locator::LocatorType
end

"""$(TYPEDSIGNATURES)
Generator function."""
function ParticleVariables(SG::SpectralGrid)
    (; architecture, nparticles, NF, ArrayType) = SG
    (; ParticleVector) = SG
    VectorNF = array_type(architecture, NF, 1)
    locator = RingGrids.AnvilLocator(NF, nparticles; architecture)
    return ParticleVariables{ParticleVector, VectorNF, typeof(locator)}(;
            nparticles, locator)
end

Adapt.@adapt_structure ParticleVariables

export DiagnosticVariables

"""All diagnostic variables.
$(TYPEDFIELDS)"""
struct DiagnosticVariables{
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
    SpectralTracerTuple,    # <: NamedTuple{Symbol, SpectralVariable3D}
    GridTracerTuple,        # <: NamedTuple{Symbol, GridVariable3D}
    ParticleVector,         # <: AbstractField
    VectorType,             # <: AbstractVector
    PhysicsTuple,           # <: NamedTyple{Symbol, ...} 
    ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
    LocatorType,            # <: AbstractLocator
    RefValueNF,             # <: Base.RefValue{NF}
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
    tendencies::Tendencies{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple}
    
    "Gridded prognostic variables"
    grid::GridVariables{GridType, GridVariable2D, GridVariable3D, GridTracerTuple}
    
    "Intermediate variables for the dynamical core"
    dynamics::DynamicsVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ScratchMemoryType}
    
    # TODO: different parameteric type here? 
    "Global fields returned from physics parameterizations, including land and ocean"
    physics::PhysicsTuple
    
    "Intermediate variables for the particle advection"
    particles::ParticleVariables{ParticleVector, VectorType, LocatorType}
    
    "Average temperature of every horizontal layer [K]"
    temp_average::VectorType

    "Scale applied to vorticity and divergence"
    scale::RefValueNF

    # full constructor that infers correct type parameters, mainly for adapt_structure / GPU etc.
    function DiagnosticVariables(
        spectrum::SpectrumType,
        grid::GridType,
        nlayers::Int,
        nparticles::Int,
        tendencies::Tendencies{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple},
        grid_variables::GridVariables{GridType, GridVariable2D, GridVariable3D, GridTracerTuple},
        dynamics::DynamicsVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ScratchMemoryType},
        physics::PhysicsTuple,
        particles::ParticleVariables{ParticleVector, VectorType, LocatorType},
        temp_average::VectorType,
        scale::RefValueNF,
    ) where {
        SpectrumType,           # <: AbstractSpectrum
        GridType,               # <: AbstractGrid
        SpectralVariable2D,     # <: LowerTriangularArray
        SpectralVariable3D,     # <: LowerTriangularArray
        GridVariable2D,         # <: AbstractField
        GridVariable3D,         # <: AbstractField
        SpectralTracerTuple,    # <: NamedTuple{Symbol, SpectralVariable3D}
        GridTracerTuple,        # <: NamedTuple{Symbol, GridVariable3D}
        ParticleVector,         # <: AbstractField
        VectorType,             # <: AbstractVector
        PhysicsTuple,
        ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
        LocatorType,            # <: AbstractLocator
        RefValueNF,             # <: Base.RefValue{NF}
    }
        return new{
            typeof(spectrum), typeof(grid), SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple,
            ParticleVector, VectorType, typeof(physics),
            typeof(dynamics.scratch_memory), typeof(particles.locator), 
            typeof(scale)
        }(
            spectrum, grid, nlayers, nparticles,
            tendencies, grid_variables, dynamics, physics, particles,
            temp_average, scale,
        )
    end

    # full constructor with parameters
    function DiagnosticVariables{
        SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple,
        ParticleVector, VectorType, PhysicsTuple,
        ScratchMemoryType, LocatorType, RefValueNF
    }(
        spectrum, grid, nlayers, nparticles,
        tendencies, grid_variables, dynamics, physics, particles,
        temp_average, scale,
    ) where {
        SpectrumType,           # <: AbstractSpectrum
        GridType,               # <: AbstractGrid
        SpectralVariable2D,     # <: LowerTriangularArray
        SpectralVariable3D,     # <: LowerTriangularArray
        GridVariable2D,         # <: AbstractField
        GridVariable3D,         # <: AbstractField
        SpectralTracerTuple,    # <: NamedTuple{Symbol, SpectralVariable3D}
        GridTracerTuple,        # <: NamedTuple{Symbol, GridVariable3D}
        ParticleVector,         # <: AbstractField
        VectorType,             # <: AbstractVector
        PhysicsTuple,
        ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
        LocatorType,            # <: AbstractLocator
        RefValueNF,             # <: Base.RefValue{NF}
    }
        return new{
            SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple,
            ParticleVector, VectorType, typeof(physics),
            ScratchMemoryType, LocatorType, RefValueNF,
        }(
            spectrum, grid, nlayers, nparticles,
            tendencies, grid_variables, dynamics, physics, particles,
            temp_average, scale,
        )
    end
end

Base.eltype(::DiagnosticVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D}) where {SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D <: AbstractField{NF}} where NF = NF 
Architectures.array_type(::DiagnosticVariables{SpectrumType, GridType, SpectralVariable2D}) where {SpectrumType, GridType, SpectralVariable2D <: LowerTriangularArray{NF, N, ArrayType}} where {NF, N, ArrayType <: AbstractArray} = nonparametric_type(ArrayType)
Architectures.array_type(::DiagnosticVariables{SpectrumType, GridType, SpectralVariable2D}) where {SpectrumType, GridType, SpectralVariable2D <: AbstractArray} = SpectralVariable2D

"""$(TYPEDSIGNATURES)
Generator function."""
function DiagnosticVariables(model::AbstractModel)

    SG = model.spectral_grid
    (; spectral_transform, tracers) = model

    if typeof(model) <: PrimitiveEquation
        nbands_shortwave = get_nbands(model.shortwave_radiation)
        nbands_longwave = get_nbands(model.longwave_radiation)
    else 
        nbands_shortwave = nbands_longwave = 0
    end

    (; spectrum, grid, nparticles, NF, nlayers) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG
    (; VectorType, ParticleVector) = SG
    nlayers_soil = model.land.nlayers 

    tendencies = Tendencies(SG; tracers)
    grid_variables = GridVariables(SG; tracers)
    dynamics = DynamicsVariables(SG; spectral_transform)
    particles = ParticleVariables(SG)
    temp_average = SG.VectorType(undef, nlayers)

    # allocate parameterization variables 
    variable_names = get_diagnostic_variables(model)

    # TODO: currently this is just a drop-in replacement, later we should have nlayers_soil from the land model and not from the SpectralGrid
    land = initialize_variables(SG, nlayers_soil, variable_names.land...)
    atmosphere = initialize_variables(SG, nlayers, variable_names.atmosphere...)
    ocean = initialize_variables(SG, 1, variable_names.ocean...)

    physics = merge(atmosphere, (ocean=ocean, land=land))

    scale = Ref(one(NF))

    SpectralTracerTuple = typeof(tendencies.tracers_tend)
    GridTracerTuple = typeof(tendencies.tracers_tend_grid)

    return DiagnosticVariables{
        typeof(spectrum), typeof(grid), SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D, SpectralTracerTuple, GridTracerTuple,
        ParticleVector, VectorType, typeof(physics),
        typeof(dynamics.scratch_memory), typeof(particles.locator), typeof(scale)
    }(
        spectrum, grid, nlayers, nparticles,
        tendencies, grid_variables, dynamics, physics, particles,
        temp_average, scale,
    )
end

Adapt.@adapt_structure DiagnosticVariables

function Base.show(
    io::IO,
    diagn::DiagnosticVariables,
)

    NF = eltype(diagn)
    ArrayType = array_type(diagn)

    println(io, "DiagnosticVariables{$NF, $ArrayType}")
    
    (; spectrum, nlayers, nparticles) = diagn
    grid = diagn.grid_used   # TODO grid is used by 'GridVariables'
    nlat = RingGrids.get_nlat(grid)
    Grid = nonparametric_type(grid)

    ntracers = length(diagn.grid.tracers_grid)
    println(io, "├ spectrum: T$(truncation(spectrum)), $nlayers layers, $ntracers tracers")
    println(io, "├ grid: $nlat-ring, $nlayers-layer $Grid")
    println(io, "├ tendencies::Tendencies")
    println(io, "├ grid::GridVariables")
    println(io, "├ dynamics::DynamicsVariables")
    println(io, "├ physics::PhysicsVariables")
    println(io, "├ particles::ParticleVariables, $nparticles particles")
    println(io, "├ temp_average::$(typeof(diagn.temp_average))")
    print(io,   "└ scale: $(diagn.scale[])")
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
