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
    lnp_vert_adv_grid   ::Grid                          # vertical advection of log surf pressure
    humid_tend_grid     ::Grid                          # specific humidity
end

function Base.zeros(::Type{Tendencies},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF
    
    @unpack Grid, nresolution = G
    @unpack lmax, mmax = S
    
    # use one more l for size compat with vector quantities
    vor_tend        = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # vorticity
    div_tend        = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # divergence
    temp_tend       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # absolute Temperature
    humid_tend      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # specific humidity
    u_tend          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # zonal velocity
    v_tend          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # meridional velocity
    u_tend_grid     = zeros(Grid{NF},nresolution)                               # zonal velocity
    v_tend_grid     = zeros(Grid{NF},nresolution)                               # meridional velocity
    temp_tend_grid  = zeros(Grid{NF},nresolution)                               # temperature
    lnp_vert_adv_grid   = zeros(Grid{NF},nresolution)                           # vert adv of log surf pres
    humid_tend_grid = zeros(Grid{NF},nresolution)                               # specific humidity

    return Tendencies(  vor_tend,div_tend,temp_tend,humid_tend,
                        u_tend,v_tend,u_tend_grid,v_tend_grid,
                        temp_tend_grid,lnp_vert_adv_grid,humid_tend_grid)
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
    U_grid              ::Grid  # zonal velocity *coslat [m/s]
    V_grid              ::Grid  # meridional velocity *coslat [m/s]
    temp_grid_anomaly   ::Grid  # absolute temperature anomaly [K]
end

function Base.zeros(::Type{GridVariables},G::Geometry{NF}) where NF

    @unpack Grid, nresolution = G
    vor_grid            = zeros(Grid{NF},nresolution)   # vorticity
    div_grid            = zeros(Grid{NF},nresolution)   # divergence
    temp_grid           = zeros(Grid{NF},nresolution)   # absolute temperature
    temp_virt_grid      = zeros(Grid{NF},nresolution)   # virtual temperature
    humid_grid          = zeros(Grid{NF},nresolution)   # specific humidity
    geopot_grid         = zeros(Grid{NF},nresolution)   # geopotential
    U_grid              = zeros(Grid{NF},nresolution)   # zonal velocity *coslat
    V_grid              = zeros(Grid{NF},nresolution)   # meridonal velocity *coslat
    temp_grid_anomaly   = zeros(Grid{NF},nresolution)   # absolute temperature anolamy

    return GridVariables(   vor_grid,div_grid,temp_grid,temp_virt_grid,humid_grid,geopot_grid,
                            U_grid,V_grid,temp_grid_anomaly)
end

"""
    DynamicsVariables{Grid<:AbstractGrid,NF<:AbstractFloat}

Struct holding intermediate quantities for the dynamics of a given layer."""
struct DynamicsVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}

    # VELOCITY VECTOR (U,V) in spectral with *coslat scaling included
    u_coslat::LowerTriangularMatrix{Complex{NF}}
    v_coslat::LowerTriangularMatrix{Complex{NF}}

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
    div_sum_above   ::Grid                          # sum of (thickness-weighted) div from top to k 
    div_weighted    ::Grid                          # = ∇⋅((uₖ,vₖ)Δpₖ), weighted by pres thick
    temp_virt       ::LowerTriangularMatrix{Complex{NF}}    # virtual temperature spectral for geopot
    geopot          ::LowerTriangularMatrix{Complex{NF}}    # geopotential on full layers

    # VERTICAL VELOCITY (̇̇dσ/dt)
    σ_tend          ::Grid                          # = dσ/dt, on half levels at k+1/2
    # σ_m             ::Grid                              # TODO what's that, also on half levels at k+1/2

    # ###------Defined in surface_pressure_tendency!()
    # u_mean             ::Array{NF,2}  # Mean gridpoint zonal velocity over all levels
    # v_mean             ::Array{NF,2}  # Mean gridpoint meridional velocity over all levels
    # div_mean           ::Array{NF,2}  # Mean gridpoint divergence over all levels

    # pres_gradient_spectral_x ::Array{Complex{NF},2} #X Gradient of the surface pressure, spectral space
    # pres_gradient_spectral_y ::Array{Complex{NF},2} #Y Gradient of the surface pressure, spectral space

    # pres_gradient_grid_x ::Array{NF,2} #X Gradient of the surface pressure, grid point space
    # pres_gradient_grid_y ::Array{NF,2} #X Gradient of the surface pressure, grid point space

    # ###------Defined in vertical_velocity_tendency!()
    # sigma_tend ::Array{NF,3} #vertical velocity in sigma coords
    # sigma_m    ::Array{NF,3} #some related quantity. What is this physically?
    # puv        ::Array{NF,3} #(ug -umean)*px + (vg -vmean)*py

    # ###------Defined in zonal_wind_tendency!()
    # sigma_u ::Array{NF,3}  #some quantity used for later calculations

    # ###------Defined in vor_div_tendency_and_corrections!()
    # L2_velocity_complex ::Array{Complex{NF},2} # -laplacian(0.5*(u**2+v**2))

    # ###-----Defined in tendencies.jl/get_spectral_tendencies!()
    # vertical_mean_divergence ::Array{Complex{NF},2}
    # sigdtc ::Array{Complex{NF},3} # what is this quantity, physically?
    # dumk ::Array{Complex{NF},3} #ditto
    # spectral_geopotential ::Array{Complex{NF},3} #This should probably go elsewhere
end

function Base.zeros(::Type{DynamicsVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF
    
    @unpack lmax, mmax = S
    @unpack Grid, nresolution = G

    # VELOCITY VECTOR (U,V) in spectral with *coslat scaling included
    u_coslat = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    v_coslat = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # MULTI-PURPOSE VECTOR (a,b), work array to be reused in various places
    a = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    b = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    a_grid = zeros(Grid{NF},nresolution)
    b_grid = zeros(Grid{NF},nresolution)

    # SHALLOW WATER and PRIMITIVE EQUATION MODEL, bernoulli = 1/2*(u^2 + v^2) + Φ
    bernoulli_grid  = zeros(Grid{NF},nresolution)
    bernoulli       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    
    # VERTICAL INTEGRATION
    uv∇lnp          = zeros(Grid{NF},nresolution)   # = (uₖ,vₖ)⋅∇ln(pₛ), pressure flux
    div_sum_above   = zeros(Grid{NF},nresolution)   # sum of thickness-weighted divs above
    div_weighted    = zeros(Grid{NF},nresolution)   # = ∇⋅((uₖ,vₖ)Δpₖ), weighted by pres thick
    temp_virt       = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    geopot          = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # VERTICAL VELOCITY (̇̇dσ/dt)
    σ_tend = zeros(Grid{NF},nresolution)    # = dσ/dt, on half levels at k+1/2
    # σ_m = zeros(Grid{NF},nresolution)       # TODO what's that?

    # sigma_tend  = zeros(NF,nlon,nlat,nlev+1)
    # sigma_m     = zeros(NF,nlon,nlat,nlev+1)
    # puv         = zeros(NF,nlon,nlat,nlev)
    # sigma_u     = zeros(NF,nlon,nlat,nlev+1)

    # L2_velocity_complex         = zeros(Complex{NF},lmax+2,mmax+1)

    # vertical_mean_divergence    = zeros(Complex{NF},lmax+2,mmax+1)
    # sigdtc                      = zeros(Complex{NF},lmax+2,mmax+1,nlev+1)
    # dumk                        = zeros(Complex{NF},lmax+2,mmax+1,nlev+1)
    # spectral_geopotential       = zeros(Complex{NF},lmax+2,mmax+1,nlev)

    return DynamicsVariables(   u_coslat,v_coslat,
                                a,b,a_grid,b_grid,
                                bernoulli_grid,bernoulli,
                                uv∇lnp,div_sum_above,div_weighted,temp_virt,geopot,
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

    U_mean_grid::Grid                               # vertical average of: zonal velocity *coslat
    V_mean_grid::Grid                               # meridional velocity *coslat
    div_mean_grid::Grid                             # divergence
    div_mean::LowerTriangularMatrix{Complex{NF}}    # divergence (in spectral though)
    
    precip_large_scale::Grid
    precip_convection::Grid

    npoints::Int    # number of grid points
end

function Base.zeros(::Type{SurfaceVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack Grid, nresolution, npoints = G
    @unpack lmax, mmax = S

    # log of surface pressure and tendency thereof
    pres_grid = zeros(Grid{NF},nresolution)
    pres_tend = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    pres_tend_grid = zeros(Grid{NF},nresolution)

    # gradients of log surface pressure
    ∇lnp_x = zeros(Grid{NF},nresolution)    # zonal gradient of log surf pressure
    ∇lnp_y = zeros(Grid{NF},nresolution)    # meridional gradient of log surf pres

    # vertical averaged (weighted by σ level thickness) velocities (*coslat) and divergence
    U_mean_grid = zeros(Grid{NF},nresolution)
    V_mean_grid = zeros(Grid{NF},nresolution)
    div_mean_grid = zeros(Grid{NF},nresolution)
    div_mean = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # precipitation fields
    precip_large_scale = zeros(Grid{NF},nresolution)
    precip_convection = zeros(Grid{NF},nresolution)

    return SurfaceVariables(pres_grid,pres_tend,pres_tend_grid,
                            ∇lnp_x,∇lnp_y,
                            U_mean_grid,V_mean_grid,div_mean_grid,div_mean,
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

# LOOP OVER ALL GRID POINTS
eachgridpoint(diagn::DiagnosticVariables) = Base.OneTo(diagn.npoints)
eachgridpoint(layer::DiagnosticVariablesLayer) = Base.OneTo(layer.npoints)
eachgridpoint(surface::SurfaceVariables) = Base.OneTo(surface.npoints)