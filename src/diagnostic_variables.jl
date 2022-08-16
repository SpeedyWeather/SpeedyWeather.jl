"""
    Tendencies{NF<:AbstractFloat}

Struct holding the tendencies of the prognostic spectral variables for a given layer."""
struct Tendencies{NF<:AbstractFloat}
    vor_tend  ::LowerTriangularMatrix{Complex{NF}}      # Vorticity of horizontal wind field [1/s]
    div_tend  ::LowerTriangularMatrix{Complex{NF}}      # Divergence of horizontal wind field [1/s]
    temp_tend ::LowerTriangularMatrix{Complex{NF}}      # Absolute temperature [K]
    humid_tend::LowerTriangularMatrix{Complex{NF}}      # Specific humidity [g/kg]
    
    u_tend    ::Matrix{NF}                              # zonal velocity
    v_tend    ::Matrix{NF}                              # meridional velocity
    humid_grid_tend::Matrix{NF}                         # specific humidity
    temp_grid_tend ::Matrix{NF}                         # temperature
end

function Base.zeros(::Type{Tendencies},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF
    
    @unpack nlon, nlat = G
    @unpack lmax, mmax = S
    
    # use one more l for size compat with vector quantities
    vor_tend    = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # vorticity
    div_tend    = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # divergence
    temp_tend   = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # absolute Temperature
    humid_tend  = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)   # specific humidity
    
    # tendencies from parametrizations in grid-point space
    u_tend          = zeros(NF,nlon,nlat)                                   # zonal velocity
    v_tend          = zeros(NF,nlon,nlat)                                   # meridional velocity
    humid_grid_tend = zeros(NF,nlon,nlat)                                   # specific humidity
    temp_grid_tend  = zeros(NF,nlon,nlat)                                   # temperature

    return Tendencies(  vor_tend,div_tend,temp_tend,humid_tend,
                        u_tend,v_tend,humid_grid_tend,temp_grid_tend)
end

"""
    GridVariables{NF<:AbstractFloat}

Struct holding the prognostic spectral variables of a given layer in grid point space."""
struct GridVariables{NF<:AbstractFloat}
    vor_grid            ::Matrix{NF}    # vorticity
    div_grid            ::Matrix{NF}    # divergence
    temp_grid           ::Matrix{NF}    # absolute temperature [K]
    humid_grid          ::Matrix{NF}    # specific_humidity
    geopot_grid         ::Matrix{NF}    # geopotential
    U_grid              ::Matrix{NF}    # zonal velocity *coslat [m/s]
    V_grid              ::Matrix{NF}    # meridional velocity *coslat [m/s]
    temp_grid_anomaly   ::Matrix{NF}    # absolute temperature anomaly [K]
end

function Base.zeros(::Type{GridVariables},G::Geometry{NF}) where NF

    @unpack nlon, nlat = G

    vor_grid            = zeros(NF,nlon,nlat) # vorticity
    div_grid            = zeros(NF,nlon,nlat) # divergence
    temp_grid           = zeros(NF,nlon,nlat) # absolute Temperature
    humid_grid          = zeros(NF,nlon,nlat) # specific humidity
    geopot_grid         = zeros(NF,nlon,nlat) # geopotential
    U_grid              = zeros(NF,nlon,nlat) # zonal velocity *coslat
    V_grid              = zeros(NF,nlon,nlat) # meridonal velocity *coslat
    temp_grid_anomaly   = zeros(NF,nlon,nlat) # absolute temperature anolamy

    return GridVariables(   vor_grid,div_grid,temp_grid,humid_grid,geopot_grid,
                            U_grid,V_grid,temp_grid_anomaly)
end

"""
    DynamicsVariables{NF<:AbstractFloat}

Struct holding intermediate quantities for the dynamics of a given layer."""
struct DynamicsVariables{NF<:AbstractFloat}

    ### VORTICITY INVERSION
    u_coslat        ::LowerTriangularMatrix{Complex{NF}}    # = U = cosθ*u, zonal velocity *cos(latitude)
    v_coslat        ::LowerTriangularMatrix{Complex{NF}}    # = V = cosθ*v, meridional velocity *cos(latitude)

    # VORTICITY ADVECTION
    uω_coslat⁻¹_grid::Matrix{NF}                            # = u(ζ+f)/coslat on the grid
    vω_coslat⁻¹_grid::Matrix{NF}                            # = v(ζ+f)/coslat on the grid
    uω_coslat⁻¹     ::LowerTriangularMatrix{Complex{NF}}    # = u(ζ+f)/coslat in spectral space
    vω_coslat⁻¹     ::LowerTriangularMatrix{Complex{NF}}    # = v(ζ+f)/coslat in spectral space

    # SHALLOW WATER
    bernoulli_grid  ::Matrix{NF}                            # bernoulli potential = 1/2(u^2+v^2) + g*η
    bernoulli       ::LowerTriangularMatrix{Complex{NF}}    # spectral bernoulli potential
    uh_coslat⁻¹_grid::Matrix{NF}                            # volume flux uh/coslat on grid
    vh_coslat⁻¹_grid::Matrix{NF}                            # volume flux vh/coslat on grid
    uh_coslat⁻¹     ::LowerTriangularMatrix{Complex{NF}}    # uh/coslat in spectral
    vh_coslat⁻¹     ::LowerTriangularMatrix{Complex{NF}}    # vh/coslat in spectral

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
    @unpack nlon, nlat = G

    # BAROTROPIC VORTICITY EQUATION (vector quantities require one more degree l)
    u_coslat = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    v_coslat = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # VORTICITY ADVECTION (vector quantities require one more degree l)
    uω_coslat⁻¹_grid = zeros(NF,nlon,nlat)
    vω_coslat⁻¹_grid = zeros(NF,nlon,nlat)
    uω_coslat⁻¹      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    vω_coslat⁻¹      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # SHALLOW WATER (bernoulli is a scalar quantity of size lmax+1,mmax+1)
    bernoulli_grid   = zeros(NF,nlon,nlat)
    bernoulli        = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    uh_coslat⁻¹_grid = zeros(NF,nlon,nlat)
    vh_coslat⁻¹_grid = zeros(NF,nlon,nlat)
    uh_coslat⁻¹      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    vh_coslat⁻¹      = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    # u_mean      = zeros(NF,nlon,nlat)           # Mean gridpoint zonal velocity over all levels
    # v_mean      = zeros(NF,nlon,nlat)           # Mean gridpoint meridional velocity over all levels
    # div_mean    = zeros(NF,nlon,nlat)           # Mean gridpoint divergence over all levels

    # # one more l for recursion in meridional gradients
    # # X,Y gradient of the surface pressure in spectral space
    # pres_gradient_spectral_x = zeros(Complex{NF},lmax+2,mmax+1)
    # pres_gradient_spectral_y = zeros(Complex{NF},lmax+2,mmax+1)

    # # X,Y gradient of the surface pressure in grid space
    # pres_gradient_grid_x = zeros(NF,nlon,nlat)
    # pres_gradient_grid_y = zeros(NF,nlon,nlat)

    # sigma_tend  = zeros(NF,nlon,nlat,nlev+1)
    # sigma_m     = zeros(NF,nlon,nlat,nlev+1)
    # puv         = zeros(NF,nlon,nlat,nlev)
    # sigma_u     = zeros(NF,nlon,nlat,nlev+1)

    # L2_velocity_complex         = zeros(Complex{NF},lmax+2,mmax+1)

    # vertical_mean_divergence    = zeros(Complex{NF},lmax+2,mmax+1)
    # sigdtc                      = zeros(Complex{NF},lmax+2,mmax+1,nlev+1)
    # dumk                        = zeros(Complex{NF},lmax+2,mmax+1,nlev+1)
    # spectral_geopotential       = zeros(Complex{NF},lmax+2,mmax+1,nlev)

    return DynamicsVariables(   u_coslat, v_coslat,
                                uω_coslat⁻¹_grid,vω_coslat⁻¹_grid,
                                uω_coslat⁻¹,vω_coslat⁻¹,

                                bernoulli_grid,bernoulli,
                                uh_coslat⁻¹_grid,vh_coslat⁻¹_grid,uh_coslat⁻¹,vh_coslat⁻¹,
                                )
                                # u_mean,v_mean,div_mean,
                                # pres_gradient_spectral_x,pres_gradient_spectral_y,
                                # pres_gradient_grid_x,pres_gradient_grid_y,
                                # sigma_tend,sigma_m,puv,sigma_u,L2_velocity_complex,
                                # vertical_mean_divergence,sigdtc,dumk,spectral_geopotential)
end

struct DiagnosticVariablesLayer{NF<:AbstractFloat}
    tendencies          ::Tendencies{NF}
    grid_variables      ::GridVariables{NF}
    dynamics_variables  ::DynamicsVariables{NF}
    npoints             ::Int       # number of grid points
end

function Base.zeros(::Type{DiagnosticVariablesLayer},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack npoints = G

    tendencies = zeros(Tendencies,G,S)
    grid_variables = zeros(GridVariables,G)
    dynamics_variables = zeros(DynamicsVariables,G,S)

    return DiagnosticVariablesLayer(tendencies,grid_variables,dynamics_variables,npoints)
end

struct SurfaceVariables{NF<:AbstractFloat}
    pres_grid::Matrix{NF}
    pres_tend::LowerTriangularMatrix{Complex{NF}}

    precip_large_scale::Matrix{NF}
    precip_convection::Matrix{NF}

    npoints::Int        # number of grid points
end

function Base.zeros(::Type{SurfaceVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack npoints, nlon, nlat = G
    @unpack lmax, mmax = S

    pres_grid = zeros(NF,nlon,nlat)
    pres_tend = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

    precip_large_scale = zeros(NF,nlon,nlat)
    precip_convection = zeros(NF,nlon,nlat)

    return SurfaceVariables(pres_grid,pres_tend,
                            precip_large_scale,precip_convection,
                            npoints)
end

"""
    DiagnosticVariables{NF<:AbstractFloat}

Struct holding the diagnostic variables."""
struct DiagnosticVariables{NF<:AbstractFloat}
    layers  ::Vector{DiagnosticVariablesLayer{NF}}
    surface ::SurfaceVariables{NF}
    nlev    ::Int       # number of vertical levels
    npoints ::Int       # number of grid points
end

function Base.zeros(::Type{DiagnosticVariables},
                    G::Geometry{NF},
                    S::SpectralTransform{NF}) where NF

    @unpack nlev,npoints = G

    layers = [zeros(DiagnosticVariablesLayer,G,S) for _ in 1:nlev]
    surface = zeros(SurfaceVariables,G,S)

    return DiagnosticVariables(layers,surface,nlev,npoints)
end

DiagnosticVariables(G::Geometry{NF},S::SpectralTransform{NF}) where NF = zeros(DiagnosticVariables,G,S)

# LOOP OVER ALL GRID POINTS
eachgridpoint(diagn::DiagnosticVariables) = 1:diagn.npoints
eachgridpoint(layer::DiagnosticVariablesLayer) = 1:layer.npoints
eachgridpoint(surface::SurfaceVariables) = 1:surface.npoints