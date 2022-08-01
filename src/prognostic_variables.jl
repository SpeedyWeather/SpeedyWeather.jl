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

struct PrognosticVariables{NF<:AbstractFloat}
    # data
    layers::Vector{PrognosticVariablesLeapfrog{NF}} # each element = 1 vertical layer (incl leapfrog dim)    
    pres::Vector{LowerTriangularMatrix{NF}}         # 2-element leapfrog vec of log of surface pressure [log(hPa)]

    # dimensions
    lmax::Int
    mmax::Int
    n_leapfrog::Int
    n_layers::Int
end

# ZERO GENERATOR FUNCTIONS
# general version
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},m::Integer,n::Integer) where NF
    vor   = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    div   = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    temp  = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    humid = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    return PrognosticVariablesLayer{NF}(vor,div,temp,humid)
end

# reduce size of unneeded variables if ModelSetup is provided
function Base.zeros(::Type{PrognosticVariablesLayer{NF}},model::ModelSetup,m::Integer,n::Integer) where NF
    vor   = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    m,n = model isa BarotropicModel ? (0,0) : (m,n)         # other variables not needed for BarotropicModel
    div   = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    m,n = model isa ShallowWaterModel ? (0,0) : (m,n)       # other variables not needed for ShallowWaterModel
    temp  = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    humid = zeros(LowerTriangularMatrix{Complex{NF}},m,n)
    return PrognosticVariablesLayer{NF}(vor,div,temp,humid)
end

# create leapfrog time steps as 2-element vector of PrognosticVariablesLayer
function Base.zeros(::Type{PrognosticVariablesLeapfrog{NF}},m::Integer,n::Integer) where NF
    return PrognosticVariablesLeapfrog{NF}(
        [zeros(PrognosticVariablesLayer{NF},m,n) for _ in 1:N_LEAPFROG])
end

# also pass on model if available
function Base.zeros(::Type{PrognosticVariablesLeapfrog{NF}},model::ModelSetup,m::Integer,n::Integer) where NF
    return PrognosticVariablesLeapfrog{NF}(
        [zeros(PrognosticVariablesLayer{NF},model,m,n) for _ in 1:N_LEAPFROG])   # 2 leapfrog time steps for vor,div,temp,humid
end

# general function to initiate all prognostic variables with zeros
function Base.zeros(::Type{PrognosticVariables{NF}},m::Integer,n::Integer,k::Integer) where NF
    layers = [zeros(PrognosticVariablesLeapfrog{NF},m,n) for _ in 1:k]  # k layers
    pres = [zeros(LowerTriangularMatrix{NF},m,n) for _ in 1:N_LEAPFROG]          # 2 leapfrog time steps for pres
    return PrognosticVariables{NF}(layers,pres,m,n,N_LEAPFROG,k)
end

# pass on model to reduce size
function Base.zeros(::Type{PrognosticVariables{NF}},model::ModelSetup,m::Integer,n::Integer,k::Integer) where NF
    layers = [zeros(PrognosticVariablesLeapfrog{NF},model,m,n) for _ in 1:k]    # k layers
    m,n = model isa BarotropicModel ? (0,0) : (m,n)                             # pressure not needed for BarotropicModel
    pres = [zeros(LowerTriangularMatrix{NF},m,n) for _ in 1:N_LEAPFROG]         # 2 leapfrog time steps for pres
    return PrognosticVariables{NF}(layers,pres,m,n,N_LEAPFROG,k)
end
