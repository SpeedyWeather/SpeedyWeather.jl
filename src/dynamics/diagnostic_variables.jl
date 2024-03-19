"""Tendencies of the prognostic variables (incl. zonal/meridional velocity) in spectral and
grid-point space. Fields are
$(TYPEDFIELDS)"""
Base.@kwdef struct Tendencies{NF<:AbstractFloat, Grid<:AbstractGrid{NF}}

    "Spectral resolution, max degree of spherical harmonics"
    trunc::Int

    "GridVariable3D-point resolution parameter,
        number of latitude rings on one hemisphere, (Eq. included)"
    nlat_half::Int

    "Number of vertical layers"
    nlayers::Int
    
    "Spectral tendency of: Vorticity of horizontal wind field [1/s]"
    vor_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Spectral tendency of: Divergence of horizontal wind field [1/s]"
    div_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Spectral tendency of: Log of surface pressure [log(Pa)],
        or interface displacement [m] for ShallowWater"
    pres_tend::LTM{Complex{NF}} = zeros(LTM{Complex{NF}}, trunc+2, trunc+1)

    "Spectral tendency of: Absolute temperature [K]"
    temp_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Spectral tendency of: Specific humidity [g/kg]"
    humid_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Spectral tendency of: zonal velocity [m/s]"
    u_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Spectral tendency of: meridional velocity [m/s]"
    v_tend::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
    
    "Grid-point tendency of: zonal velocity [m/s]"
    u_tend_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "Grid-point tendency of: meridinoal velocity [m/s]"
    v_tend_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "Grid-point tendency of: log of surface pressure [log(Pa)]"
    pres_tend_grid::Grid = zeros(Grid, nlat_half)

    "Grid-point tendency of: temperature [K]"
    temp_tend_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "Grid-point tendency of: specific humidity [kg/kg]"
    humid_tend_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
end

"""Prognostic variables in grid-point space. Fields are
$TYPEDFIELDS."""
Base.@kwdef struct GridVariables{NF<:AbstractFloat, Grid<:AbstractGrid{NF}}
    "Number of latitudes on one hemisphere (Eq. incl.), resolution parameter for any grid"
    nlat_half::Int    
    
    "Number of vertical layers"
    nlayers::Int
    
    "vorticity [1/s], but scaled with radius"
    vor_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "divergence [1/s], but scaled with radius"
    div_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "Logarithm of surface pressure [log(Pa)], or interface displacement [m] for ShallowWater"
    pres_grid::Grid = zeros(Grid, nlat_half)

    "absolute temperature [K]"
    temp_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "absolute temperature of previous time step [K]"
    temp_grid_prev::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "virtual tempereature [K]"
    temp_virt_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "specific_humidity [kg/kg]"
    humid_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "zonal velocity [m/s]"
    u_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "meridional velocity [m/s]"
    v_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "zonal velocity of previous time step [m/s]"
    u_grid_prev::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
    
    "meridional velocity of previous time step [m/s]"
    v_grid_prev::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)
end

"""
Intermediate quantities for the dynamics of a given layer.
$(TYPEDFIELDS)"""
Base.@kwdef struct DynamicsVariables{NF<:AbstractFloat, Grid<:AbstractGrid{NF}}

    "Spectral resolution, max degree of spherical harmonics"
    trunc::Int

    "Number of latitudes on one hemisphere (Eq. incl.), resolution parameter for any grid"
    nlat_half::Int    
    
    "Number of vertical layers"
    nlayers::Int

    # MULTI-PURPOSE VECTOR (a, b), work array to be reused in various places, examples:
    # uω_coslat⁻¹, vω_coslat⁻¹ = a, b                        (all models)
    # uω_coslat⁻¹_grid, vω_coslat⁻¹_grid = a_grid, b_grid    (all models)
    # uh_coslat⁻¹, vh_coslat⁻¹ = a, b                        (ShallowWaterModel)
    # uh_coslat⁻¹_grid, vh_coslat⁻¹_grid = a_grid, b_grid    (ShallowWaterModel)
    # Bernoulli potential: 1/2*(u^2+v^2) + Φ = a, a_grid     (ShallowWater + PrimitiveEquation)
    
    "Multi-purpose SpectralVariable3D a"
    a::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)

    "Multi-purpose SpectralVariable3D b"
    b::SpectralVariable3D{Complex{NF}} = 
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)

    "Multi-purpose GridVariable3D a_grid"
    a_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "Multi-purpose GridVariable3D b_grid"
    b_grid::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    
    # VERTICAL INTEGRATION
    "= (uₖ, vₖ)⋅∇ln(pₛ), pressure flux"
    uv∇lnp::GridVariable3D{NF, Grid} = 
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "sum of Δσₖ-weighted uv∇lnp above"
    uv∇lnp_sum_above::GridVariable3D{NF, Grid} =
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "sum of div_weighted from top to k"
    div_sum_above::GridVariable3D{NF, Grid} =
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "vertical velocity in sigma coordinates (dσ/dt), on half levels below, at k+1/2"
    σ_tend::GridVariable3D{NF, Grid} =
        zeros(GridVariable3D{NF, Grid}, nlat_half, nlayers)

    "Spectral virtual temperature for geopotential calculation"
    temp_virt::SpectralVariable3D{Complex{NF}} =
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)

    "Geopotential on full layers"
    geopot::SpectralVariable3D{Complex{NF}} =
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nlayers)
end

"""
Diagnostic variables for the surface layer.
$(TYPEDFIELDS)"""
Base.@kwdef struct SurfaceVariables{NF<:AbstractFloat, Grid<:AbstractGrid{NF}}

    "Spectral resolution, max degree of spherical harmonics"
    trunc::Int

    "Number of latitudes on one hemisphere (Eq. incl.), resolution parameter for any grid"
    nlat_half::Int

    # PRESSURE GRADIENT
    "zonal gradient of log surf pressure"
    ∇lnp_x::Grid = zeros(Grid, nlat_half)        
    
    "meridional gradient of log surf pres"
    ∇lnp_y::Grid = zeros(Grid, nlat_half)        

    # VERTICAL AVERAGES
    "vertical average of: zonal velocity [m/s]"
    u_mean_grid::Grid = zeros(Grid, nlat_half)   
    
    "vertical average of: meridional velocity [m/s]"
    v_mean_grid::Grid = zeros(Grid, nlat_half)   
    
    "vertical average of: divergence [1/s], but scaled by radius"
    div_mean_grid::Grid = zeros(Grid, nlat_half) 
    
    "vertical average of: divergence [1/s] in spectral space, but scaled by radius"
    div_mean::LTM{Complex{NF}} = zeros(LTM{Complex{NF}}, trunc+2, trunc+1)
    
    # PHYSICS
    "large scale precipitation (for output)"
    precip_large_scale::Grid = zeros(Grid, nlat_half)

    "convective precipitation (for output)"
    precip_convection::Grid = zeros(Grid, nlat_half)

    "cloud top [m]"
    cloud_top::Grid = zeros(Grid, nlat_half)

    "Soil moisture availability [1], fraction of surface water available to evaporate"
    soil_moisture_availability::Grid = zeros(Grid, nlat_half)
end

"""Diagnostic variables for particle advection. Fields are
$(TYPEDFIELDS)"""
Base.@kwdef struct ParticleVariables{
    NF<:AbstractFloat,
    Grid<:AbstractGrid
} <: AbstractDiagnosticVariables
    "Number of particles"
    nparticles::Int

    "Number of latitudes on one hemisphere (Eq. incld.), resolution parameter of GridVariable3D"
    nlat_half::Int

    "Work array: predicted particle locations"
    locations::Vector{Particle{NF}} = zeros(Particle{NF}, nparticles)

    "Work array: velocity u"
    u::Vector{NF} = zeros(NF, nparticles)

    "Work array: velocity v"
    v::Vector{NF} = zeros(NF, nparticles)

    "Work array: velocity w = dσ/dt"
    σ_tend::Vector{NF} = zeros(NF, nparticles)

    "Interpolator to interpolate velocity fields onto particle positions"
    interpolator::AnvilInterpolator{NF,Grid} = AnvilInterpolator(NF, Grid, nlat_half, nparticles)
end


export DiagnosticVariables

"""
All diagnostic variables.
$(TYPEDFIELDS)"""
Base.@kwdef struct DiagnosticVariables{
    NF<:AbstractFloat,
    Grid<:AbstractGrid{NF},
    Model<:ModelSetup
} <: AbstractDiagnosticVariables

    trunc::Int
    nlat_half::Int
    nlayers::Int
    npoints2D::Int
    nparticles::Int

    # DYNAMICS
    tendencies::Tendencies{NF, Grid} = 
        Tendencies{NF, Grid}(; trunc, nlat_half, nlayers)
    grid_variables::GridVariables{NF, Grid} = 
        GridVariables{NF, Grid}(; nlat_half, nlayers)
    dynamics_variables::DynamicsVariables{NF, Grid} = 
        DynamicsVariables{NF, Grid}(; trunc, nlat_half, nlayers)
    surface::SurfaceVariables{NF, Grid} =
        SurfaceVariables{NF, Grid}(; trunc, nlat_half)
    temp_average::Vector{NF} = zeros(NF, nlayers)
    scale::Base.RefValue{NF} = Ref(one(NF))
    
    # PHYSICS
    columns::Vector{ColumnVariables{NF}} = 
        [ColumnVariables{NF}(; nlev=nlayers) for thread in 1:Threads.nthreads()]

    # PARTICLES
    particles::ParticleVariables{NF, Grid} =
        ParticleVariables{NF, Grid}(; nlat_half, nparticles)
end

# generator function based on a SpectralGrid
function Base.zeros(
    ::Type{DiagnosticVariables},
    SG::SpectralGrid,
    Model::Type{<:ModelSetup}
)
    (; trunc, NF, Grid, nlat_half, nlev, npoints, nparticles) = SG
    DiagnosticVariables{NF, Grid{NF}, Model}(; trunc, nlat_half, nlayers=nlev, npoints2D=npoints, nparticles)
end

DiagnosticVariables(SG::SpectralGrid, Model::Type{<:ModelSetup}) = zeros(DiagnosticVariables, SG, Model)
DiagnosticVariables(SG::SpectralGrid, model::ModelSetup) = zeros(DiagnosticVariables, SG, model_class(model))

# LOOP OVER ALL GRID POINTS (extend from RingGrids module)
RingGrids.eachgridpoint(diagn::DiagnosticVariables) = Base.OneTo(diagn.npoints)

struct DiagnosticVariablesLayer end