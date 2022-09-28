"""One vertical layer of prognostic variables represented by their spectral coefficients."""
struct PrognosticVariablesLayer{NF<:AbstractFloat}
    # all matrices are of size lmax x mmax
    vor  ::LowerTriangularMatrix{Complex{NF}}   # Vorticity of horizontal wind field [1/s]
    div  ::LowerTriangularMatrix{Complex{NF}}   # Divergence of horizontal wind field [1/s]
    temp ::LowerTriangularMatrix{Complex{NF}}   # Absolute temperature [K]
    humid::LowerTriangularMatrix{Complex{NF}}   # Specific humidity [g/kg]
end

const N_LEAPFROG = 2

"""All vertical layers of the prognostic variables and one layer for surface pressure."""
struct PrognosticVariablesLeapfrog{NF<:AbstractFloat}
    leapfrog::Vector{PrognosticVariablesLayer{NF}}  # 2-element vector for two leapfrog time steps
end

struct SurfaceLeapfrog{NF<:AbstractFloat}
    leapfrog::Vector{LowerTriangularMatrix{Complex{NF}}}     # 2-element vector for two leapfrog time steps
end

struct PrognosticVariables{NF<:AbstractFloat, M<:ModelSetup}
    # data
    layers::Vector{PrognosticVariablesLeapfrog{NF}} # each element = 1 vertical layer (incl leapfrog dim)    
    pres::SurfaceLeapfrog{NF}                       # 2-element leapfrog vec of log of surface pressure [log(hPa)]

    # dimensions
    lmax::Int
    mmax::Int
    n_leapfrog::Int
    nlev::Int
end

# ZERO GENERATOR FUNCTIONS
# general version
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},lmax::Integer,mmax::Integer) where NF
    # use one more l for size compat with vector quantities
    vor   = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    div   = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    temp  = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    humid = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    return PrognosticVariablesLayer(vor,div,temp,humid)
end

# reduce size of unneeded variables if ModelSetup is provided
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},model::ModelSetup,lmax::Integer,mmax::Integer) where NF
    # use one more l for size compat with vector quantities
    vor = has(model, :vor) ? zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1) : LowerTriangularMatrix{Complex{NF}}(undef, 0, 0)
    div = has(model, :div) ? zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1) : LowerTriangularMatrix{Complex{NF}}(undef, 0, 0)
    temp = has(model, :temp) ? zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1) : LowerTriangularMatrix{Complex{NF}}(undef, 0, 0)
    humid = has(model, :humid) ? zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1) : LowerTriangularMatrix{Complex{NF}}(undef, 0, 0)

    return PrognosticVariablesLayer(vor,div,temp,humid)
end

# create leapfrog time steps as 2-element vector of PrognosticVariablesLayer
function Base.zeros(::Type{PrognosticVariablesLeapfrog{NF}},lmax::Integer,mmax::Integer) where NF
    # use one more l for size compat with vector quantities
    return PrognosticVariablesLeapfrog([zeros(PrognosticVariablesLayer{NF},lmax,mmax) for _ in 1:N_LEAPFROG])
end

function Base.zeros(::Type{SurfaceLeapfrog{NF}},lmax::Integer,mmax::Integer) where NF
    # use one more l for size compat with vector quantities
    return SurfaceLeapfrog([zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1) for _ in 1:N_LEAPFROG])
end

# also pass on model if available
function Base.zeros(::Type{PrognosticVariablesLeapfrog{NF}},
                    model::ModelSetup,
                    lmax::Integer,
                    mmax::Integer) where NF
    # 2 leapfrog time steps for vor,div,temp,humid
    return PrognosticVariablesLeapfrog([zeros(PrognosticVariablesLayer{NF},model,lmax,mmax) for _ in 1:N_LEAPFROG])   
end

# general function to initiate all prognostic variables with zeros
function Base.zeros(::Type{PrognosticVariables{NF}},
                    lmax::Integer,
                    mmax::Integer,
                    nlev::Integer) where NF

    layers = [zeros(PrognosticVariablesLeapfrog{NF},lmax,mmax) for _ in 1:nlev]     # k layers
    pres = zeros(SurfaceLeapfrog{NF},lmax,mmax)                                     # 2 leapfrog time steps for pres
    return PrognosticVariables{NF,ModelSetup}(layers,pres,lmax,mmax,N_LEAPFROG,nlev)
end

# pass on model to reduce size
function Base.zeros(::Type{PrognosticVariables{NF}},
                    model::ModelSetup,
                    lmax::Integer,
                    mmax::Integer,
                    nlev::Integer) where NF

    layers = [zeros(PrognosticVariablesLeapfrog{NF},model,lmax,mmax) for _ in 1:nlev]   # k layers
    pres = has(model, :pres) ? zeros(SurfaceLeapfrog{NF},lmax,mmax) : zeros(SurfaceLeapfrog{NF},-2,-1)  # 2 leapfrog time steps for pres       
    return PrognosticVariables{NF,typeof(model)}(layers,pres,lmax,mmax,N_LEAPFROG,nlev)
end

has(progn::PrognosticVariables{NF,M}, var_name::Symbol) where {NF,M} = has(M, var_name)

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
    @assert has(progn, varname)

    for (progn_layer, var_layer) in zip(progn.layers, var)
        _set_var_core!(getfield(progn_layer.leapfrog[lf], varname), var_layer)
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
             var::Vector{<:AbstractGrid}, 
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
    set_vorticity!(progn::PrognosticVariables, varargs...; kwargs...)

See [`get_var`](@ref)
"""
set_vorticity!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :vor, varargs...; kwargs...)

"""
    set_divergence!(progn::PrognosticVariables, varargs...; kwargs...)

See [`get_var`](@ref)
"""
set_divergence!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :div, varargs...; kwargs...)

"""
    set_temperature!(progn::PrognosticVariables, varargs...; kwargs...)

See [`get_var`](@ref)
"""
set_temperature!(progn::PrognosticVariables, varargs...; kwargs...) = set_var!(progn, :temp, varargs...; kwargs...)

"""
    set_humidity!(progn::PrognosticVariables, varargs...; kwargs...)

See [`get_var`](@ref)
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

    _set_var_core!(progn.pres.leapfrog[lf], pressure)

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
    @assert has(progn, var_name)
    return [getfield(layer.leapfrog[lf], var_name) for layer in progn.layers]
end 

get_vorticity(progn::PrognosticVariables; kwargs...) = get_var(progn, :vor; kwargs...)
get_divergence(progn::PrognosticVariables; kwargs...) = get_var(progn, :div; kwargs...)
get_temperature(progn::PrognosticVariables; kwargs...) = get_var(progn, :temp; kwargs...)
get_humidity(progn::PrognosticVariables; kwargs...) = get_var(progn, :humid; kwargs...)
get_pressure(progn::PrognosticVariables; lf::Integer=1) = progn.pres.leapfrog[lf]



