"""
    abstract type AbstractHEALPix4Grid{T} <: AbstractGrid{T} end

An `AbstractHEALPix4Grid` is a horizontal grid similar to the standard OctahedralGrid,
but the number of points in the ring closest to the Poles starts from 4 instead of 20,
and the longitude of the first point in each ring is shifted as in HEALPixGrid.
Also, different latitudes can be used."""
abstract type AbstractHEALPix4Grid{T} <: AbstractGrid{T} end

get_nresolution(G::AbstractHEALPix4Grid) = G.nside
get_nlat_half(::Type{<:AbstractHEALPix4Grid},nside::Integer) = nside
nlat_odd(::Type{<:AbstractHEALPix4Grid}) = true
get_nlon_max(::Type{<:AbstractHEALPix4Grid},nside::Integer) = 4nside

function get_nlon_per_ring(::Type{<:AbstractHEALPix4Grid},nside::Integer,j::Integer)
    nlat = 2nside-1
    @assert 0 < j <= nlat "Ring $j is outside J$nside grid."
    nlat_half = (nlat+1)÷2
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_healpix4(nside,j)
end

get_npoints(::Type{<:AbstractHEALPix4Grid},nside::Integer) = npoints_healpix4(nside)
get_lon(::Type{<:AbstractHEALPix4Grid},nside::Integer) = Float64[]    # only defined for full grids

function get_colatlons(Grid::Type{<:AbstractHEALPix4Grid},nside::Integer)
    npoints = get_npoints(Grid,nside)
    colat = get_colat(Grid,nside)

    colats = zeros(npoints)
    lons = zeros(npoints)

    ij = 1
    for j in 1:2nside-1
        nlon = get_nlon_per_ring(Grid,nside,j)
        lon = collect(π/nlon:2π/nlon:2π)

        colats[ij:ij+nlon-1] .= colat[j]
        lons[ij:ij+nlon-1] .= lon

        ij += nlon
    end
    return colats, lons
end

function each_index_in_ring(::Type{<:AbstractHEALPix4Grid}, # function for HEALPix4 grids
                            j::Integer,                     # ring index north to south
                            nside::Integer)                 # resolution param

    @boundscheck 0 < j < 2nside || throw(BoundsError)       # ring index valid?
    if j <= nside                                           # northern hemisphere incl Equator
        index_1st = 2j*(j-1) + 1                            # first in-ring index i
        index_end = 2j*(j+1)                                # last in-ring index i
    else                                                    # southern hemisphere 
        n = 4nside^2                                        # total number of points
        j = 2nside - j                                      # count ring index from south pole
        index_1st = n - 2j*(j+1) + 1                        # count backwards
        index_end = n - 2j*(j-1)
    end
    return index_1st:index_end                              # range of i's in ring
end


"""
    J = HEALPix4Grid{T}

A HEALPix4 grid with 4 base faces, each `nside`x`nside` grid points, each covering the same area.
The values of all grid points are stored in a vector field `v` that unravels the data 0 to 360˚,
then ring by ring, which are sorted north to south."""
struct HEALPix4Grid{T} <: AbstractHEALPix4Grid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nside::Int      # nside^2 is the number of pixel in each of the 12 base pixel

    HEALPix4Grid{T}(data,nside) where T = length(data) == npoints_healpix4(nside) ?
    new(data,nside) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create an J$nside HEALPix4Grid{$T}.")
end

npoints_healpix4(nside::Integer) = 4nside^2
nside_healpix4(npoints::Integer) = round(Int,sqrt(npoints/4))  # inverse of npoints_healpix4
nlat_healpix4(nside::Integer) = 2nside-1
nlon_healpix4(nside::Integer,j::Integer) = 4j

# infer nside from data vector length, infer parametric type from eltype of data
HEALPix4Grid{T}(data::AbstractVector) where T = HEALPix4Grid{T}(data,nside_healpix4(length(data)))
HEALPix4Grid(data::AbstractVector,n::Integer...) = HEALPix4Grid{eltype(data)}(data,n...)

truncation_order(::Type{<:HEALPix4Grid}) = 3                 # cubic (in longitude)
get_truncation(::Type{<:HEALPix4Grid},nside::Integer) = nside-1
get_resolution(::Type{<:HEALPix4Grid},trunc::Integer) = roundup_fft(trunc+1)
get_colat(G::Type{<:HEALPix4Grid},nside::Integer) =
    ( colat_half=[acos(1-(j/nside)^2) for j in 1:nside]; vcat(colat_half,π.-colat_half[end-1:-1:1]) )
# get_lon() is already implemented for AbstractHEALPix4Grid
# get_colatlons() is already implemented for AbstractHEALPix4Grid
# each_index_in_ring() is already implemented for AbstractHEALPix4Grid
full_grid(::Type{<:HEALPix4Grid}) = FullHEALPix4Grid    # the full grid with same latitudes

matrix_size(G::HEALPix4Grid) = (2G.nside,2G.nside)
matrix_size(::Type{HEALPix4Grid},nside::Integer)=(2nside,2nside)
Base.Matrix(G::HEALPix4Grid{T};kwargs...) where T = Matrix!(zeros(T,matrix_size(G)...),G;kwargs...)
"""
    Matrix!(M::AbstractMatrix,
            G::HEALPix4Grid;
            quadrant_rotation=(0,1,2,3),
            matrix_quadrant=((2,2),(1,2),(1,1),(2,1)),
            )

Sorts the gridpoints in `G` into the matrix `M` without interpolation.
Every quadrant of the grid `G` is rotated as specified in `quadrant_rotation`,
0 is no rotation, 1 is 90˚ clockwise, 2 is 180˚ etc. Grid quadrants are counted
eastward starting from 0˚E. The grid quadrants are moved into the matrix quadrant
(i,j) as specified. Defaults are equivalent to centered at 0˚E and a rotation
such that the North Pole is at M's midpoint."""
Matrix!(M::AbstractMatrix,G::HEALPix4Grid;kwargs...) = Matrix!((M,G);kwargs...)

"""
    Matrix!(MGs::Tuple{AbstractMatrix{T},HEALPix4Grid}...;kwargs...)

Like `Matrix!(::AbstractMatrix,::HEALPix4Grid)` but for simultaneous
processing of tuples `((M1,G1),(M2,G2),...)` with matrices `Mi` and grids `Gi`.
All matrices and grids have to be of the same size respectively."""
function Matrix!(   MGs::Tuple{AbstractMatrix{T},HEALPix4Grid}...;
                    quadrant_rotation::NTuple{4,Integer}=(0,1,2,3),     # = 0˚, 90˚, 180˚, 270˚ anti-clockwise
                    matrix_quadrant::NTuple{4,Tuple{Integer,Integer}}=((2,2),(1,2),(1,1),(2,1)),
                    ) where T
                    
    ntuples = length(MGs)

    # check that the first (matrix,grid) tuple has corresponding sizes
    M,G = MGs[1]
    m,n = size(M)
    @boundscheck m == n || throw(BoundsError)
    @boundscheck m == 2*G.nside || throw(BoundsError)

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
    nlat_half = G.nside         # number of latitude rings on one hemisphere incl Equator
    nside = G.nside           # side length of a basepixel matrix (=nside in HEALPix)

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
            r = min(j,nlat_half) - iq             # row in matrix m (1-based)
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
    return Tuple(Mi for (Mi,Gi) in MGs)
end