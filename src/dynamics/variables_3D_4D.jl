## GRID 3D
export GridVariable3D

"""
A horizontal grid::AbstractGrid with a third dimension (time or vertical).
As the grid is a AbstractVector (with a single index only) a GridVariable3D
is a subtype of AbstractMatrix (and not AbstractArray{T,3}) with the horizontal
dimension unravelled into its first dimension. To be generated and indexed like

    nlat_half = 24      # resolution parameter of the horizontal grid
    nlayer = 8          # number of vertical layers (or timesteps)
    G = zeros(GridVariable3D, OctahedralGaussianGrid{Float32}, nlat_half, nlayer)
    G[50, 3]            # the 50th grid point on the 3rd layer

Fields are
$(TYPEDFIELDS)"""
struct GridVariable3D{NF, Grid<:AbstractGrid{NF}} <: AbstractMatrix{NF}
    data::Vector{Grid}
    npoints2D::Int
end

function Base.zeros(
    ::Type{GridVariable3D},
    Grid::Type{<:AbstractGrid{NF}},
    nlat_half::Integer,
    k::Integer,
) where NF
    npoints2D = RingGrids.get_npoints(Grid, nlat_half)
    G3D = [zeros(Grid, nlat_half) for _ in 1:k]
    GridVariable3D{NF, Grid}(G3D, npoints2D)
end

Base.size(G3D::GridVariable3D) = (G3D.npoints2D, length(G3D.data))

@inline function Base.getindex(G::GridVariable3D, ij::Integer, k::Integer)
    @boundscheck 0 < k <= length(G.data) || throw(BoundsError(G, (ij, k)))
    @boundscheck 0 < ij <= G.npoints2D || throw(BoundsError(G, (ij, k)))
    @inbounds r = G.data[k][ij]
end

Base.@propagate_inbounds function Base.setindex!(G::GridVariable3D, x, ij::Integer, k::Integer)
    @boundscheck 0 < k <= length(G.data) || throw(BoundsError(G, (ij, k)))
    @boundscheck 0 < ij <= G.npoints2D || throw(BoundsError(G, (ij, k)))
    @inbounds G.data[k][ij] = x
end

## GRID 4D
export GridVariable4D

"""
A horizontal grid::AbstractGrid with a third and forth dimension (vertical and time).
As the grid is a AbstractVector (with a single index only) a GridVariable3D
is a subtype of AbstractArray{T, 3} (and not 4) with the horizontal dimension
unravelled into its first dimension. To be generated and indexed like

    nlat_half = 24      # resolution parameter of the horizontal grid
    nlayer = 8          # number of vertical layers
    ntimesteps = 2      # number of time steps
    G = zeros(GridVariable4D, OctahedralGaussianGrid{Float32}, nlat_half, nlayer, ntimesteps)
    G[50, 3, 1]         # the 50th grid point on the 3rd layer and the 1st timestep

Fields are
$(TYPEDFIELDS)"""
struct GridVariable4D{NF, Grid<:AbstractGrid{NF}} <: AbstractArray{NF, 3}
    data::Matrix{Grid}
    npoints2D::Int
end

function Base.zeros(
    ::Type{GridVariable4D},
    Grid::Type{<:AbstractGrid{NF}},
    nlat_half::Integer,
    k::Integer,
    l::Integer,
) where NF
    npoints2D = RingGrids.get_npoints(Grid, nlat_half)
    G4D = [zeros(Grid, nlat_half) for _ in 1:k, _ in 1:l]
    GridVariable4D{NF, Grid}(G4D, npoints2D)
end

Base.size(G4D::GridVariable4D) = (G4D.npoints2D, size(G4D.data)...)

@inline function Base.getindex(G::GridVariable4D, ij::Integer, k::Integer, l::Integer)
    @boundscheck 0 < k <= size(G.data,1) || throw(BoundsError(G, (ij, k, l)))
    @boundscheck 0 < l <= size(G.data,2) || throw(BoundsError(G, (ij, k, l)))
    @boundscheck 0 < ij <= G.npoints2D || throw(BoundsError(G, (ij, k, l)))
    @inbounds G.data[k, l][ij]
end

Base.@propagate_inbounds function Base.setindex!(G::GridVariable4D, x, ij::Integer, k::Integer, l::Integer)
    @boundscheck 0 < k <= size(G.data,1) || throw(BoundsError(G, (ij, k, l)))
    @boundscheck 0 < l <= size(G.data,2) || throw(BoundsError(G, (ij, k, l)))
    @boundscheck 0 < ij <= G.npoints2D || throw(BoundsError(G, (ij, k, l)))
    @inbounds G.data[k, l][ij] = x
end

## SPECTRAL 3D
export SpectralVariable3D

"""
A horizontal LowerTriangularMatrix for spherical harmonic coefficients with a third
dimension (time or vertical). Uses either (1) two indices lm, k with lm being the unravelled
single index for the LowerTriangularMatrix (skipping the zeros in the upper triangle)
and k the vertical layer or the timestep. Or (2) three indices l, m, k with l, m
for the (1-indexed) spherical harmonic of degree l and order m.
To be generated and indexed like

    m, n = 4, 4         # size of LowerTriangularMatrix
    nlayer = 8          # number of vertical layers (or timesteps)
    S = zeros(SpectralVariable3D{Float32}, m, n, nlayer)
    S[1, 1, 3]          # l=m=1 (1-indexed) harmonic on layer 3
    S[2, 3]             # single-indexed harmonic l=2, m=1 (1-indexed) harmonic on layer 3

Fields are
$(TYPEDFIELDS)"""
struct SpectralVariable3D{NF} <: AbstractArray{NF, 3}
    data::Vector{LowerTriangularMatrix{NF}}
    m::Int      # degrees in the LowerTriangularMatrix (1st dimension)
    n::Int      # orders in the LowerTriangularMatrix  (2nd dimension)
end

function Base.zeros(
    ::Type{SpectralVariable3D{NF}},
    m::Integer,     # length of 1st dim in LowerTriangularMatrix (degrees)
    n::Integer,     # length of 2nd dim in LowerTriangularMatrix (orders)
    k::Integer,     # length of 3rd dim (vertical or time)
) where NF
    S3D = [zeros(LowerTriangularMatrix{NF}, m, n) for _ in 1:k]
    SpectralVariable3D{NF}(S3D, m, n)
end

Base.size(S3D::SpectralVariable3D) = (G4D.m, G4D.n, length(G4D.data))

## SPECTRAL 4D

export SpectralVariable4D

"""
A horizontal LowerTriangularMatrix for spherical harmonic coefficients with a third and
forth dimension (vertical and time). Uses either (1) three indices lm, k, t with lm being the unravelled
single index for the LowerTriangularMatrix (skipping the zeros in the upper triangle)
and k the vertical layer and t the timestep. Or (2) four indices l, m, k, t with l, m
for the (1-indexed) spherical harmonic of degree l and order m.
To be generated and indexed like

    m, n = 4, 4         # size of LowerTriangularMatrix
    nlayer = 8          # number of vertical layers (or timesteps)
    ntimesteps = 2      # number of timesteps
    S = zeros(SpectralVariable3D{Float32}, m, n, nlayer, ntimesteps)
    S[1, 1, 3, 1]       # l=m=1 (1-indexed) harmonic on layer 3 at timestep 1
    S[2, 3, 1]          # single-indexed harmonic l=2, m=1 (1-indexed) harmonic on layer 3, timestep 1

Fields are
$(TYPEDFIELDS)"""
struct SpectralVariable4D{NF} <: AbstractArray{NF, 4}
    data::Matrix{LowerTriangularMatrix{NF}}
    m::Int      # degrees in the LowerTriangularMatrix (1st dimension)
    n::Int      # orders in the LowerTriangularMatrix  (2nd dimension)
end

function Base.zeros(
    ::Type{SpectralVariable3D{NF}},
    m::Integer,     # length of 1st dim in LowerTriangularMatrix (degrees)
    n::Integer,     # length of 2nd dim in LowerTriangularMatrix (orders)
    k::Integer,     # length of 3rd dim (vertical)
    t::Integer,     # length of 4th dim (time)
) where NF
    S4D = [zeros(LowerTriangularMatrix{NF}, m, n) for _ in 1:k, _ in 1:t]
    SpectralVariable4D{NF}(S4D, m, n)
end

Base.size(S4D::SpectralVariable4D) = (G4D.m, G4D.n, size(G4D.data)...)
