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

struct PrognosticVariables{NF<:AbstractFloat}
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
    vor   = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    lmax,mmax = model isa BarotropicModel   ? (-2,-1) : (lmax,mmax)   # other variables not needed for BarotropicModel
    div   = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    lmax,mmax = model isa ShallowWaterModel ? (-2,-1) : (lmax,mmax)   # other variables not needed for ShallowWaterModel
    temp  = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    humid = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
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
    return PrognosticVariables(layers,pres,lmax,mmax,N_LEAPFROG,nlev)
end

# pass on model to reduce size
function Base.zeros(::Type{PrognosticVariables{NF}},
                    model::ModelSetup,
                    lmax::Integer,
                    mmax::Integer,
                    nlev::Integer) where NF

    layers = [zeros(PrognosticVariablesLeapfrog{NF},model,lmax,mmax) for _ in 1:nlev]   # k layers
    _lmax,_mmax = model isa BarotropicModel ? (-2,-1) : (lmax,mmax)   # pressure not needed for BarotropicModel
    pres = zeros(SurfaceLeapfrog{NF},_lmax,_mmax)                     # 2 leapfrog time steps for pres
    return PrognosticVariables(layers,pres,lmax,mmax,N_LEAPFROG,nlev)
end
