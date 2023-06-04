import Base: copy!

# how many time steps have to be stored for the time integration? Leapfrog = 2 
const N_STEPS = 2

struct PrognosticVariablesLayer{NF<:AbstractFloat}
    # all matrices are of size lmax x mmax
    vor  ::LowerTriangularMatrix{Complex{NF}}   # Vorticity of horizontal wind field [1/s]
    div  ::LowerTriangularMatrix{Complex{NF}}   # Divergence of horizontal wind field [1/s]
    temp ::LowerTriangularMatrix{Complex{NF}}   # Absolute temperature [K]
    humid::LowerTriangularMatrix{Complex{NF}}   # Specific humidity [g/kg]
end

struct PrognosticVariablesSurface{NF<:AbstractFloat}
    pres::LowerTriangularMatrix{Complex{NF}}    # log of surface pressure [log(Pa)]
end

struct PrognosticLayerTimesteps{NF<:AbstractFloat}
    timesteps::Vector{PrognosticVariablesLayer{NF}}     # N_STEPS-element vector for time steps
end

struct PrognosticSurfaceTimesteps{NF<:AbstractFloat}
    timesteps::Vector{PrognosticVariablesSurface{NF}}   # N_STEPS-element vector for time steps
end

struct PrognosticVariables{NF<:AbstractFloat,M<:ModelSetup}
    # data
    layers::Vector{PrognosticLayerTimesteps{NF}}     # vector of vertical layers  
    surface::PrognosticSurfaceTimesteps{NF}

    # dimensions
    trunc::Int              # max degree of spherical harmonics
    n_steps::Int            # N_STEPS time steps that are stored
    nlev::Int               # number of vertical levels

    # scaling
    scale::Base.RefValue{NF}

    # # scaling TODO think about including more flexible scaling factors here
    # scaling_vor::Base.RefValue{NF}
    # scaling_div::Base.RefValue{NF}
    # scaling_temp::Base.RefValue{NF}
    # scaling_humid::Base.RefValue{NF}
    # scaling_pres::Base.RefValue{NF}
end

# ZERO GENERATOR FUNCTIONS
# general version
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},trunc::Integer) where NF
    # use one more l for size compatibility with vector quantities
    # size trunc+2 x trunc+1 corresponds to lmax+1 x mmax
    vor   = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    div   = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    temp  = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    humid = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    return PrognosticVariablesLayer(vor,div,temp,humid)
end

# reduce size of unneeded variables if ModelSetup is provided
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},Model::Type{<:ModelSetup},trunc::Integer) where NF
    # use one more l for size compatibility with vector quantities
    LTM = LowerTriangularMatrix{Complex{NF}}
    vor = has(Model, :vor) ? zeros(LTM,trunc+2,trunc+1) : LTM(undef, 0, 0)
    div = has(Model, :div) ? zeros(LTM,trunc+2,trunc+1) : LTM(undef, 0, 0)
    temp = has(Model, :temp) ? zeros(LTM,trunc+2,trunc+1) : LTM(undef, 0, 0)
    humid = has(Model, :humid) ? zeros(LTM,trunc+2,trunc+1) : LTM(undef, 0, 0)

    return PrognosticVariablesLayer(vor,div,temp,humid)
end

function Base.zeros(::Type{PrognosticVariablesSurface{NF}},trunc::Integer) where NF
    # use one more l for size compatibility with vector quantities
    pres = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    return PrognosticVariablesSurface(pres)
end

function Base.zeros(::Type{PrognosticVariablesSurface{NF}},Model::Type{<:ModelSetup},trunc::Integer) where NF
    # use one more l for size compatibility with vector quantities
    LTM = LowerTriangularMatrix{Complex{NF}}
    pres = has(Model, :pres) ? zeros(LTM,trunc+2,trunc+1) : LTM(undef, 0, 0)
    return PrognosticVariablesSurface(pres)
end

# create time steps as N_STEPS-element vector of PrognosticVariablesLayer
function Base.zeros(::Type{PrognosticLayerTimesteps{NF}},trunc::Integer) where NF
    return PrognosticLayerTimesteps([zeros(PrognosticVariablesLayer{NF},trunc) for _ in 1:N_STEPS])
end

function Base.zeros(::Type{PrognosticSurfaceTimesteps{NF}},trunc::Integer) where NF
    return PrognosticSurfaceTimesteps([zeros(PrognosticVariablesSurface{NF},trunc) for _ in 1:N_STEPS])
end

# also pass on model if available
function Base.zeros(::Type{PrognosticLayerTimesteps{NF}},Model::Type{<:ModelSetup},trunc::Integer) where NF
    return PrognosticLayerTimesteps([zeros(PrognosticVariablesLayer{NF},Model,trunc) for _ in 1:N_STEPS])   
end

function Base.zeros(::Type{PrognosticSurfaceTimesteps{NF}},Model::Type{<:ModelSetup},trunc::Integer) where NF
    return PrognosticSurfaceTimesteps([zeros(PrognosticVariablesSurface{NF},Model,trunc) for _ in 1:N_STEPS])   
end

# general function to initiate all prognostic variables with zeros
function Base.zeros(::Type{PrognosticVariables{NF}},trunc::Integer,nlev::Integer) where NF
    layers = [zeros(PrognosticLayerTimesteps{NF},trunc) for _ in 1:nlev]    # vector of nlev layers
    surface = zeros(PrognosticSurfaceTimesteps{NF},trunc)
    
    # initialize with scale=1, wrapped in RefValue for mutability
    scale = Ref(one(NF)) 
    return PrognosticVariables{NF,ModelSetup}(layers,surface,trunc,N_STEPS,nlev,scale)
end

# pass on model to reduce size
function Base.zeros(::Type{PrognosticVariables{NF}},Model::Type{<:ModelSetup},trunc::Integer,nlev::Integer) where NF
    
    layers = [zeros(PrognosticLayerTimesteps{NF},Model,trunc) for _ in 1:nlev]      # vector of nlev layers
    PST = PrognosticSurfaceTimesteps{NF}
    surface = has(Model, :pres) ? zeros(PST,trunc) : zeros(PST,-1)
    
    # initialize with scale=1, wrapped in RefValue for mutability
    scale = Ref(one(NF))                                                    
    return PrognosticVariables{NF,Model}(layers,surface,trunc,N_STEPS,nlev,scale)
end

has(progn::PrognosticVariables{NF,M}, var_name::Symbol) where {NF,M} = has(M, var_name)

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
    lmax,mmax = size(var_old) .- (2,1)
    var_new_trunc = spectral_truncation!(var_new, lmax+1, mmax)
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

    var_sph = [spectral(var_layer) for var_layer in var]

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

    var_grid = [spectral(var_layer, Grid) for var_layer in var]

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
function set_pressure!(progn::PrognosticVariables{NF}, 
                       pressure::LowerTriangularMatrix;
                       lf::Integer=1) where NF

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
set_pressure!(progn::PrognosticVariables, pressure::AbstractGrid, M::ModelSetup; lf::Integer=1) = set_pressure!(progn, spectral(pressure, M.spectral_transform); lf=lf)

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::AbstractGrid, 
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in grid space at leapfrog index `lf`.
"""
set_pressure!(progn::PrognosticVariables, pressure::AbstractGrid, lf::Integer=1) = set_pressure!(progn, spectral(pressure); lf=lf)

"""
    set_pressure!(progn::PrognosticVariables{NF}, 
                  pressure::AbstractMatrix, 
                  Grid::Type{<:AbstractGrid}, 
                  lf::Integer=1) where NF

Sets the prognostic variable with the surface pressure in grid space at leapfrog index `lf`.
"""
set_pressure!(progn::PrognosticVariables, pressure::AbstractMatrix, Grid::Type{<:AbstractGrid}, lf::Integer=1) = set_pressure!(progn, spectral(pressure, Grid); lf=lf)
  
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