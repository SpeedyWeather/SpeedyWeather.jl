const N_STEPS = 2   # how many time steps have to be stored for the time integration? Leapfrog = 2

function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

struct PrognosticLayerTimesteps end
struct PrognosticSurfaceTimesteps end
struct PrognosticVariablesLayer end

export PrognosticVariablesOcean
Base.@kwdef struct PrognosticVariablesOcean{
    NF<:AbstractFloat,
    Grid<:AbstractGrid{NF},
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # OCEAN VARIABLES
    "Sea surface temperature [K]"
    sea_surface_temperature::Grid = zeros(Grid, nlat_half)

    "Sea ice concentration [1]"
    sea_ice_concentration::Grid = zeros(Grid, nlat_half)
end

export PrognosticVariablesLand
Base.@kwdef struct PrognosticVariablesLand{
    NF<:AbstractFloat,
    Grid<:AbstractGrid{NF},
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # LAND VARIABLES
    "Land surface temperature [K]"
    land_surface_temperature::Grid = zeros(Grid, nlat_half)

    "Snow depth [m]"
    snow_depth::Grid = zeros(Grid, nlat_half)

    "Soil moisture layer 1, volume fraction [1]"
    soil_moisture_layer1::Grid = zeros(Grid, nlat_half)

    "Soil moisture layer 2, volume fraction [1]"
    soil_moisture_layer2::Grid = zeros(Grid, nlat_half)
end

export PrognosticVariables
Base.@kwdef struct PrognosticVariables{
    NF<:AbstractFloat,
    Grid<:AbstractGrid{NF},
    M<:ModelSetup
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "max degree of spherical harmonics (0-based)"
    trunc::Int

    "Number of latitude rings on one hemisphere (Equator excl.), resolution parameter of grids"
    nlat_half::Int

    "number of vertical layers"
    nlayers::Int

    "Number time steps that are stored, 2 for leapfrog"
    nsteps::Int

    "Number of particles for particle advection"
    nparticles::Int

    # LAYERED VARIABLES
    "Vorticity of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    vor::SpectralVariable4D{Complex{NF}} =
        zeros(SpectralVariable4D{Complex{NF}}, trunc+2, trunc+1, nlayers, nsteps)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::SpectralVariable4D{Complex{NF}} =
        zeros(SpectralVariable4D{Complex{NF}}, trunc+2, trunc+1, nlayers, nsteps)

    "Absolute temperature [K]"
    temp::SpectralVariable4D{Complex{NF}} =
        zeros(SpectralVariable4D{Complex{NF}}, trunc+2, trunc+1, nlayers, nsteps)

    "Specific humidity [kg/kg]"
    humid::SpectralVariable4D{Complex{NF}} =
        zeros(SpectralVariable4D{Complex{NF}}, trunc+2, trunc+1, nlayers, nsteps)

    # SURFACE VARIABLES
    "log of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::SpectralVariable3D{Complex{NF}} =
        zeros(SpectralVariable3D{Complex{NF}}, trunc+2, trunc+1, nsteps)

    "Ocean variables, sea surface temperature and sea ice concentration"
    ocean::PrognosticVariablesOcean{NF, Grid} = PrognosticVariablesOcean{NF, Grid}(;nlat_half)
    
    "Land variables, land surface temperature, snow and soil moisture"
    land::PrognosticVariablesLand{NF, Grid} = PrognosticVariablesLand{NF, Grid}(;nlat_half)

    "Particles for particle advection"
    particles::Vector{Particle{NF}} = zeros(Particle{NF}, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::Base.RefValue{NF} = Ref(one(NF))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::Clock = Clock()
end

function PrognosticVariables(SG::SpectralGrid, model::ModelSetup)
    (; trunc, nlat_half, nlev, Grid, NF, nparticles) = SG
    Model = model_class(model)  # strip away the parameters
    return PrognosticVariables{NF, Grid{NF}, Model}(;
            trunc, nlat_half, nlayers=nlev, nsteps=N_STEPS, nparticles)
end

function Base.show(io::IO, P::PrognosticVariables)
    ζ = P.vor.data[end,1]       # create a view on surface relative vorticity
    ζ_grid = gridded(ζ)         # to grid space
    print(io, plot(ζ_grid, title="Surface relative vorticity"))
end

has(::PrognosticVariables{NF, Grid, M}, var_name::Symbol) where {NF, Grid, M} = has(M, var_name)

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

function _set_var_core!(var_old::LowerTriangularMatrix{T}, var_new::LowerTriangularMatrix{R}) where {T, R}
    lmax, mmax = size(var_old) .- (1, 1)
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
    var_sph = [spectral(var_layer, one_more_degree=true) for var_layer in var]
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