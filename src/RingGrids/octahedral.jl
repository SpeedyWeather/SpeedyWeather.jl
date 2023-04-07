"""
    abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

An `AbstractOctahedralGrid` is a horizontal grid with 16+4i longitude
points on the latitude ring i starting with i=1 around the pole.
Different latitudes can be used, Gaussian latitudes, equi-angle latitdes, or others."""
abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

get_nlon_max(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = nlon_octahedral(nlat_half)

function get_nlon_per_ring(Grid::Type{<:AbstractOctahedralGrid},nlat_half::Integer,j::Integer)
    nlat = get_nlat(Grid,nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_octahedral(j)
end

function get_colatlons(Grid::Type{<:AbstractOctahedralGrid},nlat_half::Integer)
    
    colat = get_colat(Grid,nlat_half)
    nlat = get_nlat(Grid,nlat_half)
    
    npoints = get_npoints(Grid,nlat_half)
    colats = zeros(npoints)                 # preallocate arrays
    lons = zeros(npoints)

    ij = 1                                  # continuous index
    for j in 1:nlat                         # populate arrays ring by ring
        nlon = get_nlon_per_ring(Grid,nlat_half,j)
        lon = collect(0:2π/nlon:2π-π/nlon)

        colats[ij:ij+nlon-1] .= colat[j]
        lons[ij:ij+nlon-1] .= lon

        ij += nlon
    end

    return colats, lons
end

function each_index_in_ring(Grid::Type{<:AbstractOctahedralGrid},
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    nlat = get_nlat(Grid,nlat_half)
    @boundscheck 0 < j <= nlat || throw(BoundsError)        # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j+7) - 15                           # first in-ring index i
        index_end = 2j*(j+9)                                # last in-ring index i
    else                                                    # southern hemisphere excl Equator
        j = nlat - j + 1                                    # mirror ring index around Equator
        n = get_npoints(Grid,nlat_half) + 1                 # number of grid points + 1
        index_1st = n - 2j*(j+9)                            # count backwards
        index_end = n - (2j*(j+7) - 15)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:AbstractOctahedralGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid,nlat_half) || throw(BoundsError)

    index_end = 0
    @inbounds for j in 1:nlat_half                  # North incl Eq only
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 16 + 4j                        # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    @inbounds for (j,j_mir) in zip( nlat_half+1:nlat,       # South only
                                    nlat-nlat_half:-1:1)    # reverse index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 16 + 4j_mir                    # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end


"""
    G = OctahedralGaussianGrid{T}

An Octahedral Gaussian grid that uses `nlat` Gaussian latitudes, but a decreasing number of longitude
points per latitude ring towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20,24,28,32,...nlon-4,nlon,nlon,nlon-4,...,32,28,24,20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralGaussianGrid{T} <: AbstractOctahedralGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    # check that `nlat_half` match the vector `v` length
    OctahedralGaussianGrid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == npoints_octahedral(nlat_half,false) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralGaussianGrid{$T}.")
end

# number of points and longitudes per ring on the octahedral grid
npoints_octahedral(nlat_half::Integer,nlat_oddp::Bool) =
    nlat_oddp ? max(0,4nlat_half^2 + 32nlat_half - 16) : 4nlat_half^2 + 36nlat_half # max(0,...) needed to avoid negative array size when nlat_half==0
nlat_half_octahedral(npoints::Integer,nlat_oddp::Bool) =
    nlat_oddp ? round(Int,-4+sqrt(20 + npoints/4)) : round(Int,-9/2+sqrt((9/2)^2 + npoints/4))  # inverse
nlon_octahedral(j::Integer) = 16+4j

# infer nside from data vector length, infer parametric type from eltype of data
OctahedralGaussianGrid{T}(data::AbstractVector) where T = OctahedralGaussianGrid{T}(data,nlat_half_octahedral(length(data),false))
OctahedralGaussianGrid(data::AbstractVector,n::Integer...) = OctahedralGaussianGrid{eltype(data)}(data,n...)

nlat_odd(::Type{<:OctahedralGaussianGrid}) = false
get_npoints(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,false)
get_colat(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = get_colat(FullGaussianGrid,nlat_half)
get_quadrature_weights(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
full_grid(::Type{<:OctahedralGaussianGrid}) = FullGaussianGrid    # the full grid with same latitudes
matrix_size(G::OctahedralGaussianGrid) = (2*(4+G.nlat_half),2*(4+G.nlat_half+1))

"""
    G = OctahedralClenshawGrid{T}

An Octahedral Clenshaw grid that uses `nlat` equi-spaced latitudes. Like FullClenshawGrid, the central
latitude ring is on the Equator. Like OctahedralGaussianGrid, the number of longitude points per
latitude ring decreases towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20,24,28,32,...nlon-4,nlon,nlon,nlon-4,...,32,28,24,20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralClenshawGrid{T} <: AbstractOctahedralGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    # check that `nlat_half` match the vector `v` length
    OctahedralClenshawGrid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == npoints_octahedral(nlat_half,true) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralClenshawGrid{$T}.")
end

# infer nlat_half from data vector length, infer parametric type from eltype of data
OctahedralClenshawGrid{T}(data::AbstractVector) where T = OctahedralClenshawGrid{T}(data,
                                                        nlat_half_octahedral(length(data),true))
OctahedralClenshawGrid(data::AbstractVector,n::Integer...) = OctahedralClenshawGrid{eltype(data)}(data,n...)

nlat_odd(::Type{<:OctahedralClenshawGrid}) = true
get_npoints(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,true)
get_colat(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = get_colat(FullClenshawGrid,nlat_half)
get_quadrature_weights(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
full_grid(::Type{<:OctahedralClenshawGrid}) = FullClenshawGrid    # the full grid with same latitudes

matrix_size(G::OctahedralClenshawGrid) = (2*(4+G.nlat_half),2*(4+G.nlat_half))
matrix_size(::Type{OctahedralClenshawGrid},nlat_half::Integer) = (2*(4+nlat_half),2*(4+nlat_half))
Base.Matrix(G::OctahedralClenshawGrid{T};kwargs...) where T = Matrix!(zeros(T,matrix_size(G)...),G;kwargs...)

"""
    Matrix!(M::AbstractMatrix,
            G::OctahedralClenshawGrid;
            quadrant_rotation=(0,1,2,3),
            matrix_quadrant=((2,2),(1,2),(1,1),(2,1)),
            )

Sorts the gridpoints in `G` into the matrix `M` without interpolation.
Every quadrant of the grid `G` is rotated as specified in `quadrant_rotation`,
0 is no rotation, 1 is 90˚ clockwise, 2 is 180˚ etc. Grid quadrants are counted
eastward starting from 0˚E. The grid quadrants are moved into the matrix quadrant
(i,j) as specified. Defaults are equivalent to centered at 0˚E and a rotation
such that the North Pole is at M's midpoint."""
Matrix!(M::AbstractMatrix,G::OctahedralClenshawGrid;kwargs...) = Matrix!((M,G);kwargs...)

"""
    Matrix!(MGs::Tuple{AbstractMatrix{T},OctahedralClenshawGrid}...;kwargs...)

Like `Matrix!(::AbstractMatrix,::OctahedralClenshawGrid)` but for simultaneous
processing of tuples `((M1,G1),(M2,G2),...)` with matrices `Mi` and grids `Gi`.
All matrices and grids have to be of the same size respectively."""
function Matrix!(   MGs::Tuple{AbstractMatrix{T},OctahedralClenshawGrid}...;
                    quadrant_rotation::NTuple{4,Integer}=(0,1,2,3),     # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
                    matrix_quadrant::NTuple{4,Tuple{Integer,Integer}}=((2,2),(1,2),(1,1),(2,1)),
                    ) where T
                    
    ntuples = length(MGs)

    # check that the first (matrix,grid) tuple has corresponding sizes
    M,G = MGs[1]
    m,n = size(M)
    @boundscheck m == n || throw(BoundsError)
    @boundscheck m == 2*(4+G.nlat_half) || throw(BoundsError)

    for MG in MGs   # check that all matrices and all grids are of same size
        Mi,Gi = MG
        @boundscheck size(Mi) == size(M) || throw(BoundsError)
        @boundscheck size(Gi) == size(G) || throw(BoundsError)
    end

    for q in matrix_quadrant    # check always in 2x2
        sr,sc = q
        @boundscheck ((sr in (1,2)) && (sc in (1,2))) || throw(BoundsError)
    end

    rings = eachring(G)         # index ranges for all rings
    nlat_half = G.nlat_half     # number of latitude rings on one hemisphere incl Equator
    nside = 4+G.nlat_half       # side length of a basepixel matrix

    # sort grid indices from G into matrix M
    # 1) loop over each grid point per ring
    # 2) determine quadrant (0,1,2,3) via modulo
    # 3) get longitude index iq within quadrant
    # 4) determine corresponding indices r,c in matrix M

    @inbounds for (j,ring) in enumerate(rings)
        nlon = length(ring)                         # number of grid points in ring
        for ij in ring                              # continuous index in grid
            i = ij-ring[1]                          # 0-based index in ring
            grid_quadrant = floor(Int,mod(4*i/nlon,4))  # either 0,1,2,3
            iq = i - grid_quadrant*(nlon÷4)             # 0-based index i relative to quadrant
            r = 4+min(j,nlat_half) - iq             # row in matrix m (1-based)
            c = (iq+1) + max(0,j-nlat_half)         # column in matrix m (1-based)

            # rotate indices in quadrant
            r,c = rotate_matrix_indices(r,c,nside,quadrant_rotation[grid_quadrant+1])

            # shift grid quadrant to matrix quadrant
            sr,sc = matrix_quadrant[grid_quadrant+1]
            r += (sr-1)*nside                       # shift row into matrix quadrant
            c += (sc-1)*nside                       # shift column into matrix quadrant

            for (Mi,Gi) in MGs                      # for every (matrix,grid) tuple
                Mi[r,c] = convert(T,Gi[ij])         # convert data and copy over
            end
        end
    end

    ntuples == 1 && return M
    return Tuple(Mi for (Mi,Gi) in MGs)             # 
end