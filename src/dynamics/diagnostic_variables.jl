# const LTM = LowerTriangularMatrix     # already defined in prognostic_variables

"""
Tendencies of the prognostic spectral variables for a given layer.
$(TYPEDFIELDS)"""
Base.@kwdef struct Tendencies{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    nlat_half::Int
    trunc::Int
    vor_tend  ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # Vorticity of horizontal wind field [1/s]
    div_tend  ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # Divergence of horizontal wind field [1/s]
    temp_tend ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # Absolute temperature [K]
    humid_tend::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # Specific humidity [g/kg]
    u_tend    ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # zonal velocity (spectral)
    v_tend    ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   # meridional velocity (spectral)
    u_tend_grid     ::Grid = zeros(Grid,nlat_half)                      # zonal velocity (grid)
    v_tend_grid     ::Grid = zeros(Grid,nlat_half)                      # meridinoal velocity (grid)
    temp_tend_grid  ::Grid = zeros(Grid,nlat_half)                      # temperature
    humid_tend_grid ::Grid = zeros(Grid,nlat_half)                      # specific humidity
end

# generator function based on a SpectralGrid
function Base.zeros(::Type{Tendencies},
                    SG::SpectralGrid)
    (;NF, trunc, Grid, nlat_half) = SG
    return Tendencies{NF,Grid{NF}}(;nlat_half,trunc)
end

"""
Transformed prognostic variables (plus a few others) into grid-point space.
$TYPEDFIELDS."""
Base.@kwdef struct GridVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    nlat_half::Int                                      # resolution parameter for any grid
    vor_grid        ::Grid = zeros(Grid,nlat_half)  # vorticity
    div_grid        ::Grid = zeros(Grid,nlat_half)  # divergence
    temp_grid       ::Grid = zeros(Grid,nlat_half)  # absolute temperature [K]
    temp_grid_prev  ::Grid = zeros(Grid,nlat_half)  # absolute temperature of previous time step [K]
    temp_virt_grid  ::Grid = zeros(Grid,nlat_half)  # virtual tempereature [K]  
    humid_grid      ::Grid = zeros(Grid,nlat_half)  # specific_humidity
    geopot_grid     ::Grid = zeros(Grid,nlat_half)  # geopotential (is that needed?)
    u_grid          ::Grid = zeros(Grid,nlat_half)  # zonal velocity *coslat [m/s]
    v_grid          ::Grid = zeros(Grid,nlat_half)  # meridional velocity *coslat [m/s]
    u_grid_prev     ::Grid = zeros(Grid,nlat_half)  # zonal velocity *coslat of previous time step [m/s]
    v_grid_prev     ::Grid = zeros(Grid,nlat_half)  # meridional velocity *coslat of previous time step [m/s]
end

# generator function based on a SpectralGrid
function Base.zeros(::Type{GridVariables},SG::SpectralGrid)
    (;NF, Grid, nlat_half) = SG
    return GridVariables{NF,Grid{NF}}(;nlat_half)
end

"""
Intermediate quantities for the dynamics of a given layer.
$(TYPEDFIELDS)"""
Base.@kwdef struct DynamicsVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}

    nlat_half::Int
    trunc::Int

    # MULTI-PURPOSE VECTOR (a,b), work array to be reused in various places, examples:
    # uω_coslat⁻¹, vω_coslat⁻¹ = a,b                        (all models)
    # uω_coslat⁻¹_grid, vω_coslat⁻¹_grid = a_grid,b_grid    (all models)
    # uh_coslat⁻¹, vh_coslat⁻¹ = a,b                        (ShallowWaterModel)
    # uh_coslat⁻¹_grid, vh_coslat⁻¹_grid = a_grid, b_grid   (ShallowWaterModel)
    # Bernoulli potential: 1/2*(u^2+v^2) + Φ = a,a_grid     (ShallowWater + PrimitiveEquation)
    a::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
    b::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
    a_grid::Grid = zeros(Grid,nlat_half)
    b_grid::Grid = zeros(Grid,nlat_half)
    
    # VERTICAL INTEGRATION
    uv∇lnp          ::Grid = zeros(Grid,nlat_half)  # = (uₖ,vₖ)⋅∇ln(pₛ), pressure flux
    uv∇lnp_sum_above::Grid = zeros(Grid,nlat_half)  # sum of Δσₖ-weighted uv∇lnp above
    div_sum_above   ::Grid = zeros(Grid,nlat_half)  # sum of div_weighted from top to k

    # virtual temperature spectral for geopot, geopotential on full layers
    temp_virt   ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
    geopot      ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)

    # VERTICAL VELOCITY (̇̇dσ/dt)
    σ_tend::Grid = zeros(Grid,nlat_half)            # = dσ/dt, on half levels below, at k+1/2
end

# generator function based on a SpectralGrid
function Base.zeros(::Type{DynamicsVariables},
                    SG::SpectralGrid)
    (;NF, trunc, Grid, nlat_half) = SG
    return DynamicsVariables{NF,Grid{NF}}(;nlat_half,trunc)
end

"""
All diagnostic variables for a given layer: tendencies, prognostic varibles on the grid,
and intermediate dynamics variables.
$(TYPEDFIELDS)"""
struct DiagnosticVariablesLayer{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    npoints             ::Int                   # number of grid points
    k                   ::Int                   # which vertical model level?
    
    tendencies          ::Tendencies{NF,Grid}
    grid_variables      ::GridVariables{NF,Grid}
    dynamics_variables  ::DynamicsVariables{NF,Grid}
    temp_average        ::Base.RefValue{NF}     # average temperature for this level
end

# generator function based on a SpectralGrid
function Base.zeros(::Type{DiagnosticVariablesLayer},
                    SG::SpectralGrid,
                    k::Integer=0)                   # use k=0 (i.e. unspecified) as default
    (;npoints) = SG
    tendencies = zeros(Tendencies,SG)
    grid_variables = zeros(GridVariables,SG)
    dynamics_variables = zeros(DynamicsVariables,SG)
    temp_average = Ref(zero(SG.NF))
    return DiagnosticVariablesLayer(npoints,k,tendencies,grid_variables,dynamics_variables,temp_average)
end

"""
Diagnostic variables for the surface layer.
$(TYPEDFIELDS)"""
Base.@kwdef struct SurfaceVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF}}

    nlat_half::Int
    trunc::Int
    npoints::Int

    # log surface pressure, tendency of it and gridded tendency
    pres_grid::Grid = zeros(Grid,nlat_half)
    pres_tend::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
    pres_tend_grid::Grid = zeros(Grid,nlat_half)

    ∇lnp_x::Grid = zeros(Grid,nlat_half)        # zonal gradient of log surf pressure
    ∇lnp_y::Grid = zeros(Grid,nlat_half)        # meridional gradient of log surf pres

    u_mean_grid::Grid = zeros(Grid,nlat_half)   # vertical average of: zonal velocity *coslat
    v_mean_grid::Grid = zeros(Grid,nlat_half)   # meridional velocity *coslat
    div_mean_grid::Grid = zeros(Grid,nlat_half) # divergence
    div_mean::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)    # divergence (in spectral though)
    
    precip_large_scale::Grid = zeros(Grid,nlat_half)    # large scale precipitation (for output)
    precip_convection::Grid = zeros(Grid,nlat_half)     # convective precipitation (for output)
    soil_moisture_availability::Grid = zeros(Grid,nlat_half)
end

# generator function based on a SpectralGrid
function Base.zeros(::Type{SurfaceVariables},
                    SG::SpectralGrid)

    (;NF, trunc, Grid, nlat_half, npoints) = SG
    return SurfaceVariables{NF,Grid{NF}}(;nlat_half,trunc,npoints)
end

"""
All diagnostic variables.
$(TYPEDFIELDS)"""
struct DiagnosticVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF},Model<:ModelSetup}
    layers  ::Vector{DiagnosticVariablesLayer{NF,Grid}}
    surface ::SurfaceVariables{NF,Grid}
    columns ::Vector{ColumnVariables{NF}}

    nlat_half::Int      # resolution parameter of any Grid
    nlev    ::Int       # number of vertical levels
    npoints ::Int       # number of grid points

    scale::Base.RefValue{NF}   # vorticity and divergence are scaled by radius
end

# generator function based on a SpectralGrid
function Base.zeros(
    ::Type{DiagnosticVariables},
    SG::SpectralGrid,
    Model::Type{<:ModelSetup}
)

    (;NF,Grid,nlat_half, nlev, npoints) = SG
    layers = [zeros(DiagnosticVariablesLayer,SG,k) for k in 1:nlev]
    surface = zeros(SurfaceVariables,SG)
    
    # create one column variable per thread to avoid race conditions
    nthreads = Threads.nthreads()
    columns = [ColumnVariables{NF}(;nlev) for _ in 1:nthreads]

    scale = Ref(convert(SG.NF,SG.radius))

    return DiagnosticVariables{NF,Grid{NF},Model}(
        layers,surface,columns,nlat_half,nlev,npoints,scale)
end

DiagnosticVariables(SG::SpectralGrid) = zeros(DiagnosticVariables,SG,DEFAULT_MODEL)
DiagnosticVariables(SG::SpectralGrid,Model::Type{<:ModelSetup}) = zeros(DiagnosticVariables,SG,Model)
DiagnosticVariables(SG::SpectralGrid,model::ModelSetup) = zeros(DiagnosticVariables,SG,model_class(model))

# LOOP OVER ALL GRID POINTS (extend from RingGrids module)
RingGrids.eachgridpoint(diagn::DiagnosticVariables) = Base.OneTo(diagn.npoints)
RingGrids.eachgridpoint(layer::DiagnosticVariablesLayer) = Base.OneTo(layer.npoints)
RingGrids.eachgridpoint(surface::SurfaceVariables) = Base.OneTo(surface.npoints)