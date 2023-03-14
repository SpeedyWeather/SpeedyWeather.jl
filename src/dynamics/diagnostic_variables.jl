"""
    Tendencies{Grid<:AbstractGrid,NF<:AbstractFloat}

Struct holding the tendencies of the prognostic spectral variables for a given layer."""
struct Tendencies{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    vor_tend  ::LowerTriangularMatrix{Complex{NF}}      # Vorticity of horizontal wind field [1/s]
    div_tend  ::LowerTriangularMatrix{Complex{NF}}      # Divergence of horizontal wind field [1/s]
    temp_tend ::LowerTriangularMatrix{Complex{NF}}      # Absolute temperature [K]
    humid_tend::LowerTriangularMatrix{Complex{NF}}      # Specific humidity [g/kg]
    u_tend    ::LowerTriangularMatrix{Complex{NF}}      # zonal velocity (spectral)
    v_tend    ::LowerTriangularMatrix{Complex{NF}}      # meridional velocity (spectral)
    u_tend_grid         ::Grid                          # zonal velocity (grid)
    v_tend_grid         ::Grid                          # meridinoal velocity (grid)
    temp_tend_grid      ::Grid                          # temperature
    humid_tend_grid     ::Grid                          # specific humidity
end

function Base.zeros(::Type{Tendencies},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF
    
    @unpack Grid, nlat_half = G
    @unpack lmax, mmax = S
    
    # use one more l for size compat with vector quantities
    vor_tend        = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # vorticity
    div_tend        = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # divergence
    temp_tend       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # absolute Temperature
    humid_tend      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # specific humidity
    u_tend          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # zonal velocity
    v_tend          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # meridional velocity
    u_tend_grid     = zeros(Grid{NF},nlat_half)                               # zonal velocity
    v_tend_grid     = zeros(Grid{NF},nlat_half)                               # meridional velocity
    temp_tend_grid  = zeros(Grid{NF},nlat_half)                               # temperature
    humid_tend_grid = zeros(Grid{NF},nlat_half)                               # specific humidity

    return Tendencies(  vor_tend,div_tend,temp_tend,humid_tend,
                        u_tend,v_tend,u_tend_grid,v_tend_grid,
                        temp_tend_grid,humid_tend_grid)
end

"""
    GridVariables{NF<:AbstractFloat}

Struct holding the prognostic spectral variables of a given layer in grid point space."""
struct GridVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    vor_grid            ::Grid  # vorticity
    div_grid            ::Grid  # divergence
    temp_grid           ::Grid  # absolute temperature [K]
    temp_virt_grid      ::Grid  # virtual tempereature [K]  
    humid_grid          ::Grid  # specific_humidity
    geopot_grid         ::Grid  # geopotential (is that needed?)
    u_grid              ::Grid  # zonal velocity *coslat [m/s]
    v_grid              ::Grid  # meridional velocity *coslat [m/s]
end

function Base.zeros(::Type{GridVariables},G::Geometry{NF}) where NF

    @unpack Grid, nlat_half = G
    vor_grid            = zeros(Grid{NF},nlat_half)   # vorticity
    div_grid            = zeros(Grid{NF},nlat_half)   # divergence
    temp_grid           = zeros(Grid{NF},nlat_half)   # absolute temperature
    temp_virt_grid      = zeros(Grid{NF},nlat_half)   # virtual temperature
    humid_grid          = zeros(Grid{NF},nlat_half)   # specific humidity
    geopot_grid         = zeros(Grid{NF},nlat_half)   # geopotential
    u_grid              = zeros(Grid{NF},nlat_half)   # zonal velocity *coslat
    v_grid              = zeros(Grid{NF},nlat_half)   # meridonal velocity *coslat

    return GridVariables(   vor_grid,div_grid,temp_grid,temp_virt_grid,humid_grid,geopot_grid,
                            u_grid,v_grid)
end

"""
    DynamicsVariables{Grid<:AbstractGrid,NF<:AbstractFloat}

Struct holding intermediate quantities for the dynamics of a given layer."""
struct DynamicsVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}

    # VELOCITY VECTOR (U,V) in spectral with *coslat scaling included
    U::LowerTriangularMatrix{Complex{NF}}
    V::LowerTriangularMatrix{Complex{NF}}

    # MULTI-PURPOSE VECTOR (a,b), work array to be reused in various places
    # uω_coslat⁻¹, vω_coslat⁻¹ = a,b                        (all models)
    # uω_coslat⁻¹_grid, vω_coslat⁻¹_grid = a_grid,b_grid    (all models)
    # uh_coslat⁻¹, vh_coslat⁻¹ = a,b                        (ShallowWaterModel)
    # uh_coslat⁻¹_grid, vh_coslat⁻¹_grid = a_grid, b_grid   (ShallowWaterModel)
    a::LowerTriangularMatrix{Complex{NF}}
    b::LowerTriangularMatrix{Complex{NF}}
    a_grid::Grid
    b_grid::Grid
    
    # SHALLOW WATER and PRIMITIVE EQUATION MODEL
    bernoulli_grid  ::Grid                                  # bernoulli potential = 1/2(u^2+v^2) + g*η
    bernoulli       ::LowerTriangularMatrix{Complex{NF}}    # spectral bernoulli potential

    # VERTICAL INTEGRATION
    uv∇lnp          ::Grid                          # = (uₖ,vₖ)⋅∇ln(pₛ), pressure flux
    uv∇lnp_sum_above::Grid                          # sum of Δσₖ-weighted uv∇lnp above
    div_sum_above   ::Grid                          # sum of div_weighted from top to k
    temp_virt       ::LowerTriangularMatrix{Complex{NF}}    # virtual temperature spectral for geopot
    geopot          ::LowerTriangularMatrix{Complex{NF}}    # geopotential on full layers

    # VERTICAL VELOCITY (̇̇dσ/dt)
    σ_tend          ::Grid                          # = dσ/dt, on half levels below, at k+1/2
end

function Base.zeros(::Type{DynamicsVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF
    
    @unpack lmax, mmax = S
    @unpack Grid, nlat_half = G

    # VELOCITY VECTOR (U,V) in spectral with *coslat scaling included
    U = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    V = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # MULTI-PURPOSE VECTOR (a,b), work array to be reused in various places
    a = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    b = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    a_grid = zeros(Grid{NF},nlat_half)
    b_grid = zeros(Grid{NF},nlat_half)

    # SHALLOW WATER and PRIMITIVE EQUATION MODEL, bernoulli = 1/2*(u^2 + v^2) + Φ
    bernoulli_grid  = zeros(Grid{NF},nlat_half)
    bernoulli       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    
    # VERTICAL INTEGRATION
    uv∇lnp          = zeros(Grid{NF},nlat_half)   # = (uₖ,vₖ)⋅∇ln(pₛ), pressure flux
    uv∇lnp_sum_above= zeros(Grid{NF},nlat_half)   # sum of Δσₖ-weighted uv∇lnp from 1:k-1
    div_sum_above   = zeros(Grid{NF},nlat_half)   # sum of Δσₖ-weighted div from 1:k-1
    temp_virt       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    geopot          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # VERTICAL VELOCITY (̇̇dσ/dt)
    σ_tend = zeros(Grid{NF},nlat_half)    # = dσ/dt, on half levels below, at k+1/2

    return DynamicsVariables(   U,V,
                                a,b,a_grid,b_grid,
                                bernoulli_grid,bernoulli,
                                uv∇lnp,uv∇lnp_sum_above,div_sum_above,
                                temp_virt,geopot,
                                σ_tend,
                                )
end

struct DiagnosticVariablesLayer{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    tendencies          ::Tendencies{NF,Grid}
    grid_variables      ::GridVariables{NF,Grid}
    dynamics_variables  ::DynamicsVariables{NF,Grid}
    npoints             ::Int       # number of grid points
    k                   ::Int       # which vertical model level?
end

function Base.zeros(::Type{DiagnosticVariablesLayer},
                    G::Geometry{NF},
                    S::SpectralTransform{NF};
                    k::Integer=0) where NF      # use k=0 (i.e. unspecified) as default

    @unpack npoints = G

    tendencies = zeros(Tendencies,G,S)
    grid_variables = zeros(GridVariables,G)
    dynamics_variables = zeros(DynamicsVariables,G,S)

    return DiagnosticVariablesLayer(tendencies,grid_variables,dynamics_variables,npoints,k)
end

struct SurfaceVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    pres_grid::Grid                                 # log surface pressure
    pres_tend::LowerTriangularMatrix{Complex{NF}}   # tendency of it
    pres_tend_grid::Grid                            # gridded tendency

    ∇lnp_x::Grid                                    # zonal gradient of log surf pressure
    ∇lnp_y::Grid                                    # meridional gradient of log surf pres

    u_mean_grid::Grid                               # vertical average of: zonal velocity *coslat
    v_mean_grid::Grid                               # meridional velocity *coslat
    div_mean_grid::Grid                             # divergence
    div_mean::LowerTriangularMatrix{Complex{NF}}    # divergence (in spectral though)
    
    precip_large_scale::Grid
    precip_convection::Grid

    npoints::Int    # number of grid points
end

function Base.zeros(::Type{SurfaceVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack Grid, nlat_half, npoints = G
    @unpack lmax, mmax = S

    # log of surface pressure and tendency thereof
    pres_grid = zeros(Grid{NF},nlat_half)
    pres_tend = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    pres_tend_grid = zeros(Grid{NF},nlat_half)

    # gradients of log surface pressure
    ∇lnp_x = zeros(Grid{NF},nlat_half)    # zonal gradient of log surf pressure
    ∇lnp_y = zeros(Grid{NF},nlat_half)    # meridional gradient of log surf pres

    # vertical averaged (weighted by σ level thickness) velocities (*coslat) and divergence
    u_mean_grid = zeros(Grid{NF},nlat_half)
    v_mean_grid = zeros(Grid{NF},nlat_half)
    div_mean_grid = zeros(Grid{NF},nlat_half)
    div_mean = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # precipitation fields
    precip_large_scale = zeros(Grid{NF},nlat_half)
    precip_convection = zeros(Grid{NF},nlat_half)

    return SurfaceVariables(pres_grid,pres_tend,pres_tend_grid,
                            ∇lnp_x,∇lnp_y,
                            u_mean_grid,v_mean_grid,div_mean_grid,div_mean,
                            precip_large_scale,precip_convection,
                            npoints)
end

"""
    DiagnosticVariables{Grid<:AbstractGrid,NF<:AbstractFloat}

Struct holding the diagnostic variables."""
struct DiagnosticVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    layers  ::Vector{DiagnosticVariablesLayer{NF,Grid}}
    surface ::SurfaceVariables{NF,Grid}
    nlev    ::Int       # number of vertical levels
    npoints ::Int       # number of grid points
end

function Base.zeros(::Type{DiagnosticVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack nlev,npoints = G

    layers = [zeros(DiagnosticVariablesLayer,G,S;k) for k in 1:nlev]
    surface = zeros(SurfaceVariables,G,S)

    return DiagnosticVariables(layers,surface,nlev,npoints)
end

DiagnosticVariables(G::Geometry{NF},S::SpectralTransform{NF}) where NF = zeros(DiagnosticVariables,G,S)

# LOOP OVER ALL GRID POINTS (extend from RingGrids module)
RingGrids.eachgridpoint(diagn::DiagnosticVariables) = Base.OneTo(diagn.npoints)
RingGrids.eachgridpoint(layer::DiagnosticVariablesLayer) = Base.OneTo(layer.npoints)
RingGrids.eachgridpoint(surface::SurfaceVariables) = Base.OneTo(surface.npoints)