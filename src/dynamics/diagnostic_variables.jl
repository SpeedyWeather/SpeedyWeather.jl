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
    (; spectrum, grid, nlayers) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return Tendencies{
        typeof(spectrum), typeof(grid),
        SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D,
    }(;
        spectrum, grid, nlayers
    )
end

function Adapt.adapt_structure(to, tendencies::Tendencies)
    # empty dictionaries don't get adapted per default...
    if isempty(tendencies.tracers_tend)
        vor_tend_adapted = adapt_structure(to, tendencies.vor_tend)
        u_tend_grid_adapted = adapt_structure(to, tendencies.u_tend_grid)

        return Tendencies(
            adapt_structure(to, tendencies.spectrum),
            adapt_structure(to, tendencies.grid),
            tendencies.nlayers,
            vor_tend_adapted,
            adapt_structure(to, tendencies.div_tend),
            adapt_structure(to, tendencies.temp_tend),
            adapt_structure(to, tendencies.humid_tend),
            adapt_structure(to, tendencies.u_tend),
            adapt_structure(to, tendencies.v_tend),
            adapt_structure(to, tendencies.pres_tend),
            Dict{Symbol, typeof(vor_tend_adapted)}(),
            u_tend_grid_adapted,
            adapt_structure(to, tendencies.v_tend_grid),
            adapt_structure(to, tendencies.temp_tend_grid),
            adapt_structure(to, tendencies.humid_tend_grid),
            adapt_structure(to, tendencies.pres_tend_grid),
            Dict{Symbol, typeof(u_tend_grid_adapted)}(),
        )
    else
        return Tendencies(
            adapt_structure(to, tendencies.spectrum),
            adapt_structure(to, tendencies.grid),
            tendencies.nlayers,
            adapt_structure(to, tendencies.vor_tend),
            adapt_structure(to, tendencies.div_tend),
            adapt_structure(to, tendencies.temp_tend),
            adapt_structure(to, tendencies.humid_tend),
            adapt_structure(to, tendencies.u_tend),
            adapt_structure(to, tendencies.v_tend),
            adapt_structure(to, tendencies.pres_tend),
            adapt_structure(to, tendencies.tracers_tend),
            adapt_structure(to, tendencies.u_tend_grid),
            adapt_structure(to, tendencies.v_tend_grid),
            adapt_structure(to, tendencies.temp_tend_grid),
            adapt_structure(to, tendencies.humid_tend_grid),
            adapt_structure(to, tendencies.pres_tend_grid),
            adapt_structure(to, tendencies.tracers_tend_grid),
        )
    end
end

export GridVariables

"""Transformed prognostic variables (and u, v, temp_virt) into grid-point space.
$TYPEDFIELDS."""
@kwdef struct GridVariables{
    GridType,               # <:AbstractGrid
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
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
    "Surface pressure [Pa] at previous time step (not logarithm!)"
    pres_grid_prev  ::GridVariable2D = zeros(GridVariable2D, grid)
    "Tracers [?] at previous time step"
    tracers_grid_prev::Dict{Symbol, GridVariable3D} = Dict{Symbol, GridVariable3D}()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function GridVariables(SG::SpectralGrid)
    (; grid, nlayers) = SG
    (; GridVariable2D, GridVariable3D) = SG

    return GridVariables{typeof(grid), GridVariable2D, GridVariable3D}(;
            grid, nlayers,
        )
end

function Adapt.adapt_structure(to, grid_variables::GridVariables)
    # empty dictionaries don't get adapted per default...
    if isempty(grid_variables.tracers_grid)
        vor_grid_adapted = adapt_structure(to, grid_variables.vor_grid) 

        return GridVariables(
            adapt_structure(to, grid_variables.grid), 
            grid_variables.nlayers,
            vor_grid_adapted,
            adapt_structure(to, grid_variables.div_grid), 
            adapt_structure(to, grid_variables.temp_grid), 
            adapt_structure(to, grid_variables.temp_virt_grid), 
            adapt_structure(to, grid_variables.humid_grid), 
            adapt_structure(to, grid_variables.u_grid), 
            adapt_structure(to, grid_variables.v_grid), 
            adapt_structure(to, grid_variables.pres_grid), 
            Dict{Symbol, typeof(vor_grid_adapted)}(),
            adapt_structure(to, grid_variables.random_pattern), 
            adapt_structure(to, grid_variables.temp_grid_prev), 
            adapt_structure(to, grid_variables.humid_grid_prev), 
            adapt_structure(to, grid_variables.u_grid_prev), 
            adapt_structure(to, grid_variables.v_grid_prev), 
            adapt_structure(to, grid_variables.pres_grid_prev), 
            Dict{Symbol, typeof(vor_grid_adapted)}(),
        )   
    else 
        return GridVariables(
            adapt_structure(to, grid_variables.grid), 
            grid_variables.nlayers,
            adapt_structure(to, grid_variables.vor_grid), 
            adapt_structure(to, grid_variables.div_grid), 
            adapt_structure(to, grid_variables.temp_grid), 
            adapt_structure(to, grid_variables.temp_virt_grid), 
            adapt_structure(to, grid_variables.humid_grid), 
            adapt_structure(to, grid_variables.u_grid), 
            adapt_structure(to, grid_variables.v_grid), 
            adapt_structure(to, grid_variables.pres_grid), 
            adapt_structure(to, grid_variables.tracers_grid), 
            adapt_structure(to, grid_variables.random_pattern), 
            adapt_structure(to, grid_variables.temp_grid_prev), 
            adapt_structure(to, grid_variables.humid_grid_prev), 
            adapt_structure(to, grid_variables.u_grid_prev), 
            adapt_structure(to, grid_variables.pres_grid_prev), 
            adapt_structure(to, grid_variables.tracers_grid_prev),
        )
    end
end

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
    
    "Virtual temperature [K], spectral, used for geopotential"
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
    (; architecture, spectrum, grid, nlayers) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG

    if isnothing(spectral_transform)
        return DynamicsVariables{typeof(spectrum), typeof(grid),
            SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, 
            SpeedyTransforms.ScratchMemory{NF, array_type(architecture, Complex{NF}, 3), 
            array_type(architecture, NF, 1), array_type(architecture, Complex{NF}, 1)}}(;
                spectrum, grid, nlayers,
            )
    else 
        scratch_memory = spectral_transform.scratch_memory 

        return DynamicsVariables{typeof(spectrum), typeof(grid),
            SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, typeof(scratch_memory)}(;
                spectrum, grid, nlayers, scratch_memory
            )
    end 
end

Adapt.@adapt_structure DynamicsVariables

export DiagnosticVariablesOcean
@kwdef struct DiagnosticVariablesOcean{
    GridType,
    GridVariable2D,
} <: AbstractDiagnosticVariables

    "Grid used for fields"
    grid::GridType

    "Surface sensible heat flux [W/m²], positive up"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface humidity flux [kg/s/m²], positive up"
    surface_humidity_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface shortwave radiative flux up [W/m²]"
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface longwave radiative flux up [W/m²]"
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, grid)

    "Albedo over ocean (but defined everywhere) [1]"
    albedo::GridVariable2D = zeros(GridVariable2D, grid)
end

DiagnosticVariablesOcean(SG::SpectralGrid) =
    DiagnosticVariablesOcean{typeof(SG.grid), SG.GridVariable2D}(; SG.grid)

Adapt.@adapt_structure DiagnosticVariablesOcean

export DiagnosticVariablesLand
@kwdef struct DiagnosticVariablesLand{
    GridType,
    GridVariable2D,
} <: AbstractDiagnosticVariables

    "Grid used for fields"
    grid::GridType

    "Surface sensible heat flux [W/m²], positive up"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface humidity flux [W/m²], positive up"
    surface_humidity_flux::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface shortwave radiative flux up [W/m²]"
    surface_shortwave_up::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface longwave radiative flux up [W/m²]"
    surface_longwave_up::GridVariable2D = zeros(GridVariable2D, grid)

    "Albedo over land (but defined everywhere) [1]"
    albedo::GridVariable2D = zeros(GridVariable2D, grid)

    "Availability of soil moisture to evaporation (or condensation) [1]"
    soil_moisture_availability::GridVariable2D = zeros(GridVariable2D, grid)

    "River runoff [m/s], diagnostic overflow from soil moisture"
    river_runoff::GridVariable2D = zeros(GridVariable2D, grid)
end

DiagnosticVariablesLand(SG::SpectralGrid) =
    DiagnosticVariablesLand{typeof(SG.grid), SG.GridVariable2D}(; SG.grid)

Adapt.@adapt_structure DiagnosticVariablesLand

export PhysicsVariables

"""
Diagnostic variables of the physical parameterizations.
$(TYPEDFIELDS)"""
@kwdef struct PhysicsVariables{
    GridType,               # <:AbstractGrid
    GridVariable2D,         # <: AbstractField
} <: AbstractDiagnosticVariables

    grid::GridType          # resolution of grid

    ocean::DiagnosticVariablesOcean{GridType, GridVariable2D}
    land::DiagnosticVariablesLand{GridType, GridVariable2D}

    # PRECIPITATION
    "Accumulated large-scale rain [m]"
    rain_large_scale::GridVariable2D = zeros(GridVariable2D, grid)

    "Accumulated convective rain [m]"
    rain_convection::GridVariable2D = zeros(GridVariable2D, grid)

    "Accumulated large-scale snow [m]"
    snow_large_scale::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Accumulated convective snow [m]"
    snow_convection::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Rate of total precipitation (rain+snow) [kg/m²/s]"
    total_precipitation_rate::GridVariable2D = zeros(GridVariable2D, grid)

    "Cloud top [m]"
    cloud_top::GridVariable2D = zeros(GridVariable2D, grid)      
    
    # SURFACE FLUXES
    "Surface wind speed [m/s]"
    surface_wind_speed::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface air density [kg/m³]"
    surface_air_density::GridVariable2D = zeros(GridVariable2D, grid)

    "Boundary layer drag coefficient [1]"
    boundary_layer_drag::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface air temperature [K]"
    surface_air_temperature::GridVariable2D = zeros(GridVariable2D, grid)

    "Sensible heat flux [W/m²], positive up"
    sensible_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)
    
    "Surface humidity flux [kg/s/m^2], positive up"
    surface_humidity_flux::GridVariable2D = zeros(GridVariable2D, grid)

    "Surface latent heat flux [W/m²], positive up"
    surface_latent_heat_flux::GridVariable2D = zeros(GridVariable2D, grid)

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
    (; grid, GridVariable2D) = SG

    ocean = DiagnosticVariablesOcean(SG)
    land = DiagnosticVariablesLand(SG)

    return PhysicsVariables{typeof(grid), GridVariable2D}(; grid, ocean, land)
end

Adapt.@adapt_structure PhysicsVariables

export ParticleVariables

"""Diagnostic variables for the particle advection
$(TYPEDFIELDS)"""
@kwdef struct ParticleVariables{
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
    (; NF, architecture, nparticles, ParticleVector) = SG
    VectorNF = array_type(architecture, NF, 1)
    interpolator = RingGrids.AnvilInterpolator(SG.grid, nparticles; NF)
    return ParticleVariables{ParticleVector, VectorNF, typeof(interpolator)}(;
            nparticles, interpolator)
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
    ParticleVector,         # <: AbstractField
    VectorType,             # <: AbstractVector
    ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
    Interpolator,           # <:AbstractInterpolator
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
    tendencies::Tendencies{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D}
    
    "Gridded prognostic variables"
    grid::GridVariables{GridType, GridVariable2D, GridVariable3D}
    
    "Intermediate variables for the dynamical core"
    dynamics::DynamicsVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ScratchMemoryType}
    
    "Global fields returned from physics parameterizations"
    physics::PhysicsVariables{GridType, GridVariable2D}
    
    "Intermediate variables for the particle advection"
    particles::ParticleVariables{ParticleVector, VectorType, Interpolator}
    
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
        tendencies::Tendencies{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D},
        grid_variables::GridVariables{GridType, GridVariable2D, GridVariable3D},
        dynamics::DynamicsVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D, GridVariable3D, ScratchMemoryType},
        physics::PhysicsVariables{GridType, GridVariable2D},
        particles::ParticleVariables{ParticleVector, VectorType, Interpolator},
        temp_average::VectorType,
        scale::RefValueNF,
    ) where {
        SpectrumType,           # <: AbstractSpectrum
        GridType,               # <: AbstractGrid
        SpectralVariable2D,     # <: LowerTriangularArray
        SpectralVariable3D,     # <: LowerTriangularArray
        GridVariable2D,         # <: AbstractField
        GridVariable3D,         # <: AbstractField
        ParticleVector,         # <: AbstractField
        VectorType,             # <: AbstractVector
        ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
        Interpolator,           # <: AbstractInterpolator
        RefValueNF,             # <: Base.RefValue{NF}
    }
        return new{
            typeof(spectrum), typeof(grid), SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, ParticleVector, VectorType,
            typeof(dynamics.scratch_memory), typeof(particles.interpolator), 
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
        GridVariable2D, GridVariable3D, ParticleVector, VectorType,
        ScratchMemoryType, Interpolator, RefValueNF
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
        ParticleVector,         # <: AbstractField
        VectorType,             # <: AbstractVector
        ScratchMemoryType,      # <: ArrayType{Complex{NF}, 3}
        Interpolator,           # <: AbstractInterpolator
        RefValueNF,             # <: Base.RefValue{NF}
    }
        return new{
            SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D,
            GridVariable2D, GridVariable3D, ParticleVector, VectorType,
            ScratchMemoryType, Interpolator, RefValueNF,
        }(
            spectrum, grid, nlayers, nparticles,
            tendencies, grid_variables, dynamics, physics, particles,
            temp_average, scale,
        )
    end
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

Base.eltype(::DiagnosticVariables{SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D}) where {SpectrumType, GridType, SpectralVariable2D, SpectralVariable3D, GridVariable2D <: AbstractField{NF}} where NF = NF 
Architectures.array_type(::DiagnosticVariables{SpectrumType, GridType, SpectralVariable2D}) where {SpectrumType, GridType, SpectralVariable2D <: LowerTriangularArray{NF, N, ArrayType}} where {NF, N, ArrayType <: AbstractArray} = nonparametric_type(ArrayType)

"""$(TYPEDSIGNATURES)
Generator function. If a `transform` is handed over, the same scratch memory is used."""
function DiagnosticVariables(
    SG::SpectralGrid;
    spectral_transform::Union{Nothing, SpectralTransform} = nothing,
    nbands_shortwave::Integer = 0,
    nbands_longwave::Integer = 0,
)
    (; spectrum, grid, nparticles, NF, nlayers) = SG
    (; SpectralVariable2D, SpectralVariable3D) = SG
    (; GridVariable2D, GridVariable3D) = SG
    (; VectorType, ParticleVector) = SG

    tendencies = Tendencies(SG)
    grid_variables = GridVariables(SG)
    dynamics = DynamicsVariables(SG; spectral_transform)
    physics = PhysicsVariables(SG)
    particles = ParticleVariables(SG)
    temp_average = SG.VectorType(undef, nlayers)

    scale = Ref(one(NF))

    return DiagnosticVariables{
        typeof(spectrum), typeof(grid), SpectralVariable2D, SpectralVariable3D,
        GridVariable2D, GridVariable3D, ParticleVector, VectorType,
        typeof(dynamics.scratch_memory), typeof(particles.interpolator), typeof(scale)
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

function add!(diagn::DiagnosticVariables{
    SpectrumType,           # <: AbstractSpectrum
    GridType,               # <:AbstractGrid
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractField
    GridVariable3D,         # <: AbstractField
},
    tracers::Tracer...,
    ) where {
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
