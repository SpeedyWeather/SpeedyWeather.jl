abstract type AbstractVariables end

const DEFAULT_DATE = DateTime(2000,1,1)

"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes
$(TYPEDFIELDS)."""
Base.@kwdef mutable struct Clock
    "current model time"
    time::DateTime = DEFAULT_DATE
    
    "number of days to integrate for, set in run!(::Simulation)"
    n_days::Float64 = 0

    "number of time steps to integrate for, set in initialize!(::Clock,::TimeStepper)"
    n_timesteps::Int = 0  
end

# pretty printing
function Base.show(io::IO,C::Clock)
    println(io,"$(typeof(C))")
    keys = propertynames(C)
    print_fields(io,C,keys)
end

"""
$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` in the `time_stepping`."""
function initialize!(clock::Clock,time_stepping::TimeStepper)
    clock.n_timesteps = ceil(Int,24*clock.n_days/(value(Hour(time_stepping.Δt_sec))))
    return clock
end

"""
$(TYPEDSIGNATURES)
Create and initialize a clock from `time_stepping`"""
function Clock(time_stepping::TimeStepper;kwargs...)
    clock = Clock(;kwargs...)
    initialize!(clock,time_stepping)
end

# how many time steps have to be stored for the time integration? Leapfrog = 2 
const N_STEPS = 2
const LTM = LowerTriangularMatrix       # just because it's shorter here

"""A layer of the prognostic variables in spectral space.
$(TYPEDFIELDS)"""
Base.@kwdef struct PrognosticVariablesLayer{NF<:AbstractFloat} <: AbstractVariables

    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int

    "Vorticity of horizontal wind field [1/s]"
    vor  ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)   

    "Divergence of horizontal wind field [1/s]"
    div  ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)

    "Absolute temperature [K]"
    temp ::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)

    "Specific humidity [kg/kg]"
    humid::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
end

# generator function based on spectral grid
PrognosticVariablesLayer(SG::SpectralGrid) = PrognosticVariablesLayer{SG.NF}(trunc=SG.trunc)

function Base.show(io::IO,A::AbstractVariables)
    println(io,"$(typeof(A))")
    keys = propertynames(A)
    print_fields(io,A,keys)
end

"""Collect the n time steps of PrognosticVariablesLayer
of an n-step time integration (leapfrog=2) into a single struct.
$(TYPEDFIELDS).""" 
struct PrognosticLayerTimesteps{NF<:AbstractFloat} <: AbstractVariables
    timesteps::Vector{PrognosticVariablesLayer{NF}}     # N_STEPS-element vector for time steps
end

# generator function based on spectral grid
function PrognosticLayerTimesteps(SG::SpectralGrid)
    return PrognosticLayerTimesteps([PrognosticVariablesLayer(SG) for _ in 1:N_STEPS])
end

"""The spectral and gridded prognostic variables at the surface.
$(TYPEDFIELDS)"""
Base.@kwdef struct PrognosticVariablesSurface{NF<:AbstractFloat} <: AbstractVariables

    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int

    "log of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::LTM{Complex{NF}} = zeros(LTM{Complex{NF}},trunc+2,trunc+1)
end

# generator function based on a SpectralGrid
PrognosticVariablesSurface(SG::SpectralGrid) = PrognosticVariablesSurface{SG.NF}(trunc=SG.trunc)

Base.@kwdef mutable struct PrognosticVariablesOcean{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractVariables

    "Resolution parameter of grid"
    const nlat_half::Int

    "Current time of the ocean variables"
    time::DateTime = DEFAULT_DATE

    # SEA
    "Sea surface temperature [K]"
    const sea_surface_temperature::Grid = zeros(Grid,nlat_half)

    "Sea ice concentration [1]"
    const sea_ice_concentration::Grid = zeros(Grid,nlat_half)
end

# generator function based on a SpectralGrid
function PrognosticVariablesOcean(SG::SpectralGrid)
    (;nlat_half,Grid,NF) = SG
    return PrognosticVariablesOcean{NF,Grid{NF}}(;nlat_half)
end

Base.@kwdef mutable struct PrognosticVariablesLand{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractVariables

    "Resolution parameter of grid"
    const nlat_half::Int

    "Current time of the land variables"
    time::DateTime = DEFAULT_DATE

    # LAND
    "Land surface temperature [K]"
    const land_surface_temperature::Grid = zeros(Grid,nlat_half)

    "Snow depth [m]"
    const snow_depth::Grid = zeros(Grid,nlat_half)

    "Soil moisture layer 1, volume fraction [1]"
    const soil_moisture_layer1::Grid = zeros(Grid,nlat_half)

    "Soil moisture layer 2, volume fraction [1]"
    const soil_moisture_layer2::Grid = zeros(Grid,nlat_half)
end

# generator function based on a SpectralGrid
function PrognosticVariablesLand(SG::SpectralGrid)
    (;nlat_half,Grid,NF) = SG
    return PrognosticVariablesLand{NF,Grid{NF}}(;nlat_half)
end

"""Collect the n time steps of PrognosticVariablesSurface
of an n-step time integration (leapfrog=2) into a single struct.
$(TYPEDFIELDS).""" 
struct PrognosticSurfaceTimesteps{NF<:AbstractFloat} <: AbstractVariables
    timesteps::Vector{PrognosticVariablesSurface{NF}}   # N_STEPS-element vector for time steps
end

# generator function based on spectral grid
function PrognosticSurfaceTimesteps(SG::SpectralGrid)
    return PrognosticSurfaceTimesteps([PrognosticVariablesSurface(SG) for _ in 1:N_STEPS])
end

struct PrognosticVariables{NF<:AbstractFloat,Grid<:AbstractGrid{NF},M<:ModelSetup} <: AbstractVariables

    # dimensions
    trunc::Int              # max degree of spherical harmonics
    nlat_half::Int          # resolution parameter of grids
    nlev::Int               # number of vertical levels
    n_steps::Int            # N_STEPS time steps that are stored

    layers::Vector{PrognosticLayerTimesteps{NF}}     # vector of vertical layers  
    surface::PrognosticSurfaceTimesteps{NF}
    ocean::PrognosticVariablesOcean{NF,Grid}
    land::PrognosticVariablesLand{NF,Grid}

    # scaling
    scale::Base.RefValue{NF}

    clock::Clock
end

function PrognosticVariables(SG::SpectralGrid,model::ModelSetup)
    
    (;trunc,nlat_half,nlev,Grid,NF) = SG

    # data structs
    layers = [PrognosticLayerTimesteps(SG) for _ in 1:nlev]      # vector of nlev layers
    surface = PrognosticSurfaceTimesteps(SG)
    ocean = PrognosticVariablesOcean(SG)
    land = PrognosticVariablesLand(SG)

    scale = Ref(one(NF))        # initialize with scale=1, wrapped in RefValue for mutability
    clock = Clock()

    Model = model_class(model)  # strip away the parameters
    return PrognosticVariables{NF,Grid{NF},Model}(  trunc,nlat_half,nlev,N_STEPS,
                                                    layers,surface,ocean,land,scale,clock)
end

has(::PrognosticVariables{NF,Grid,M}, var_name::Symbol) where {NF,Grid,M} = has(M, var_name)

"""
    copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

Copies entries of `progn_old` into `progn_new`. Only copies those variables that are present 
in the model of both `progn_new` and `progn_old`.
"""
function Base.copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

    var_names = propertynames(progn_old.layers[1].timesteps[1])

    for var_name in var_names
        if has(progn_new, var_name) 
            var = get_var(progn_old, var_name) 
            set_var!(progn_new, var_name, var)
        end
    end 
    pres = get_pressure(progn_old)
    set_pressure!(progn_new, pres)
    
    # synchronize the clock
    progn_new.clock.time = progn_old.clock.time

    return progn_new
end

# SET_VAR FUNCTIONS TO ASSIGN NEW VALUES TO PrognosticVariables

"""
    set_var!(progn::PrognosticVariables{NF},        
             varname::Symbol, 
             var::Vector{<:LowerTriangularMatrix};
             lf::Integer=1) where NF

Sets the prognostic variable with the name `varname` in all layers at leapfrog index `lf` 
with values given in `var` a vector with all information for all layers in spectral space.
"""
function set_var!(progn::PrognosticVariables{NF}, 
                  varname::Symbol, 
                  var::Vector{<:LowerTriangularMatrix};
                  lf::Integer=1) where NF

    @assert length(var) == length(progn.layers)
    @assert has(progn, varname) "PrognosticVariables has no variable $varname"

    for (progn_layer, var_layer) in zip(progn.layers, var)
        _set_var_core!(getfield(progn_layer.timesteps[lf], varname), var_layer)
    end 

    return progn 
end 

function _set_var_core!(var_old::LowerTriangularMatrix{T}, var_new::LowerTriangularMatrix{R}) where {T,R}
    lmax,mmax = size(var_old) .- (1,1)
    var_new_trunc = spectral_truncation!(var_new, mmax+1, mmax)
    copyto!(var_old, var_new_trunc)
end 

"""
    set_var!(progn::PrognosticVariables{NF},        
             varname::Symbol, 
             var::Vector{<:AbstractGrid};
             lf::Integer=1) where NF

Sets the prognostic variable with the name `varname` in all layers at leapfrog index `lf` 
with values given in `var` a vector with all information for all layers in grid space.
"""
function set_var!(progn::PrognosticVariables{NF}, 
                  varname::Symbol, 
                  var::Vector{<:AbstractGrid};
                  lf::Integer=1) where NF

    @assert length(var) == length(progn.layers)
    var_sph = [spectral(var_layer,one_more_degree=true) for var_layer in var]
    return set_var!(progn, varname, var_sph; lf=lf)
end 

"""
    set_var!(progn::PrognosticVariables{NF}, 
             varname::Symbol, 
             var::Vector{<:AbstractGrid}, 
             M::ModelSetup;
             lf::Integer=1) where NF

Sets the prognostic variable with the name `varname` in all layers at leapfrog index `lf` 
with values given in `var` a vector with all information for all layers in grid space.
"""
function set_var!(progn::PrognosticVariables{NF}, 
                  varname::Symbol, 
                  var::Vector{<:AbstractGrid}, 
                  M::ModelSetup;
                  lf::Integer=1) where NF

    @assert length(var) == length(progn.layers)
    
    var_sph = [spectral(var_layer, M.spectral_transform) for var_layer in var]

    return set_var!(progn, varname, var_sph; lf=lf)
end 

"""
    set_var!(progn::PrognosticVariables{NF}, 
             varname::Symbol, 
             var::Vector{<:AbstractMatrix}, 
             Grid::Type{<:AbstractGrid}=FullGaussianGrid;
             lf::Integer=1) where NF

Sets the prognostic variable with the name `varname` in all layers at leapfrog index `lf` 
with values given in `var` a vector with all information for all layers in grid space.
"""
function set_var!(progn::PrognosticVariables{NF}, 
                  varname::Symbol, 
                  var::Vector{<:AbstractMatrix}, 
                  Grid::Type{<:AbstractGrid}=FullGaussianGrid;
                  lf::Integer=1) where NF

    @assert length(var) == length(progn.layers)

    var_grid = [spectral(var_layer; Grid, one_more_degree=true) for var_layer in var]

    return set_var!(progn, varname, var_grid; lf=lf)
end 

"""
    function set_var!(progn::PrognosticVariables{NF}, 
                      varname::Symbol, 
                      s::Number;
                      lf::Integer=1) where NF

Sets all values of prognostic variable `varname` at leapfrog index `lf` to the scalar `s`.
"""
function set_var!(progn::PrognosticVariables{NF}, 
                  varname::Symbol, 
                  s::Number;
                  lf::Integer=1) where NF

    for progn_layer in progn.layers
        fill!(getfield(progn_layer.timesteps[lf], varname), s)
    end 

    return progn 
end 

"""
    set_vorticity!(progn::PrognosticVariables, varargs...; kwargs...)

See [`set_var!`](@ref)
"""
set_vorticity!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :vor, varargs...; kwargs...)

"""
    set_divergence!(progn::PrognosticVariables, varargs...; kwargs...)

See [`set_var!`](@ref)
"""
set_divergence!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :div, varargs...; kwargs...)

"""
    set_temperature!(progn::PrognosticVariables, varargs...; kwargs...)

See [`set_var!`](@ref)
"""
set_temperature!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :temp, varargs...; kwargs...)

"""
    set_humidity!(progn::PrognosticVariables, varargs...; kwargs...)

See [`set_var!`](@ref)
"""
set_humidity!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :humid, varargs...; kwargs...)

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::LowerTriangularMatrix;
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in spectral space at leapfrog index `lf`.
"""
function set_pressure!(progn::PrognosticVariables,
                       pressure::LowerTriangularMatrix;
                       lf::Integer=1)

    _set_var_core!(progn.surface.timesteps[lf].pres, pressure)

    return progn
end

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::AbstractGrid, 
                  M::ModelSetup;
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in grid space at leapfrog index `lf`.
"""
set_pressure!(progn::PrognosticVariables, pressure::AbstractGrid, M::ModelSetup; lf::Integer=1) =
    set_pressure!(progn, spectral(pressure, M.spectral_transform); lf)

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::AbstractGrid, 
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in grid space at leapfrog index `lf`.
"""
set_pressure!(progn::PrognosticVariables, pressure::AbstractGrid; lf::Integer=1) =
    set_pressure!(progn, spectral(pressure, one_more_degree=true); lf)

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::AbstractMatrix, 
                  Grid::Type{<:AbstractGrid}, 
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in grid space at leapfrog index `lf`.
"""
set_pressure!(progn::PrognosticVariables, pressure::AbstractMatrix; lf::Integer=1,
    Grid::Type{<:AbstractGrid}=FullGaussianGrid) = set_pressure!(progn, spectral(pressure; Grid, one_more_degree=true); lf)
  
"""
    get_var(progn::PrognosticVariables, var_name::Symbol; lf::Integer=1)

Returns the prognostic variable `var_name` at leapfrog index `lf` as a `Vector{LowerTriangularMatrices}`.
"""
function get_var(progn::PrognosticVariables, var_name::Symbol; lf::Integer=1)
    @assert has(progn, var_name) "PrognosticVariables has no variable $var_name"
    return [getfield(layer.timesteps[lf], var_name) for layer in progn.layers]
end 

get_vorticity(progn::PrognosticVariables; kwargs...) = get_var(progn, :vor; kwargs...)
get_divergence(progn::PrognosticVariables; kwargs...) = get_var(progn, :div; kwargs...)
get_temperature(progn::PrognosticVariables; kwargs...) = get_var(progn, :temp; kwargs...)
get_humidity(progn::PrognosticVariables; kwargs...) = get_var(progn, :humid; kwargs...)
get_pressure(progn::PrognosticVariables; lf::Integer=1) = progn.surface.timesteps[lf].pres

function Base.show(io::IO, P::PrognosticVariables)
    ζ = P.layers[end].timesteps[1].vor          # create a view on surface relative vorticity
    ζ_grid = gridded(ζ)                         # to grid space
    print(io,plot(ζ_grid,title="Surface relative vorticity"))
end