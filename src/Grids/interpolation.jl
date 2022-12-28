struct GridGeometry{G<:AbstractGrid}
    # general precomputed arrays describing the grid
    nlat_half::Int                  # number of latitude rings on one hemisphere (Eq. incl)
    nlat::Int                       # total number of latitude rings
    npoints::Int                    # total number of grid points
    latd::Vector{Float64}           # latitude of each ring, incl north pole 90˚N, ..., south pole -90˚N
    lons::Vector{Float64}           # longitudes of every grid point

    rings::Vector{UnitRange{Int}}   # for every ring a i:j unit range for ring-based indices
    nlons::Vector{Int}              # number of longitudinal points per ring
    lon_offsets::Vector{Float64}    # longitude offsets of first grid point per ring
end

GridGeometry(grid::AbstractGrid) = GridGeometry(typeof(grid),grid.nlat_half)

function GridGeometry(  Grid::Type{<:AbstractGrid},
                        nlat_half::Integer)

    nlat = get_nlat(Grid,nlat_half)                 # total number of latitude rings
    npoints = get_npoints(Grid,nlat_half)           # total number of grid points

    # LATITUDES
    colat = get_colat(Grid,nlat_half)               # colatitude in radians
    lat = π/2 .- colat                              # latitude in radians
    latd = lat*360/2π                               # 90˚...-90˚, in degrees
    latd_poles = cat(90,latd,-90,dims=1)            # latd, but poles incl
    latd_poles[end] = latd_poles[end] - eps(latd_poles[end])

    # COORDINATES for every grid point in ring order
    _,londs = get_colatlons(Grid,nlat_half)         # in radians
    londs *= (360/2π)                               # in degrees [0˚...360˚E]

    # RINGS and LONGITUDE OFFSETS
    rings = eachring(Grid,nlat_half)
    nlons = get_nlons(Grid,nlat_half,both_hemispheres=true)
    lon_offsets = [londs[ring[1]] for ring in rings]

    return GridGeometry{Grid}(nlat_half,nlat,npoints,latd_poles,londs,rings,nlons,lon_offsets)
end
    
struct InterpolationLocator{NF<:AbstractFloat}
    npoints::Int            # number of points to interpolate onto

    # to the coordinates respective indices
    js::Vector{Int}         # ring indices j such that [j,j+1) contains the point
    ij_as::Vector{Int}      # pixel index ij for top left point a on ring j
    ij_bs::Vector{Int}      # pixel index ij for top right point b on ring j
    ij_cs::Vector{Int}      # pixel index ij for bottom left point c on ring j+1
    ij_ds::Vector{Int}      # pixel index ij for bottom right point d on ring j+1

    # distances to adjacent grid points
    Δys::Vector{NF}         # distance fractions between rings
    Δabs::Vector{NF}        # distance fractions between a,b
    Δcds::Vector{NF}        # distance fractions between c,d
end

InterpolationLocator(npoints::Integer) = InterpolationLocator(Float64,npoints)

function InterpolationLocator(  ::Type{NF},
                                npoints::Integer) where {NF<:AbstractFloat}

    # to the coordinates respective indices
    js = zeros(Int,npoints)     # ring indices j such that [j,j+1) contains the point
    ij_as = zeros(Int,npoints)  # pixel index ij for top left point a on ring j
    ij_bs = zeros(Int,npoints)  # pixel index ij for top right point b on ring j
    ij_cs = zeros(Int,npoints)  # pixel index ij for bottom left point c on ring j+1
    ij_ds = zeros(Int,npoints)  # pixel index ij for bottom right point d on ring j+1

    # distances to adjacent grid points
    Δabs = zeros(NF,npoints)    # distance fractions of point on ring j-1
    Δcds = zeros(NF,npoints)    # distance fractions of point on ring j
    Δys = zeros(NF,npoints)     # distance fractions of point from j-1 to j

    return InterpolationLocator{NF}(npoints,js,ij_as,ij_bs,ij_cs,ij_ds,Δabs,Δcds,Δys)
end

struct Interpolator{NF<:AbstractFloat,G<:AbstractGrid}
    geometry::GridGeometry{G}
    locator::InterpolationLocator{NF}
end

function Interpolator(::Type{NF},grid::AbstractGrid,n::Integer) where {NF<:AbstractFloat}
    geometry = GridGeometry(grid)
    locator = InterpolationLocator(NF,n)
    return Interpolator{NF,typeof(grid)}(geometry,locator)
end

Interpolator(grid::AbstractGrid,n::Integer) = Interpolator(Float64,grid,n)

## FUNCTIONS
interpolate(θ::Real,λ::Real,A::AbstractGrid) = interpolate([θ],[λ],A)

function interpolate(θs::Vector{NF},λs::Vector{NF},A::AbstractGrid{NF2}) where {NF,NF2}
    n = length(θs)
    @assert n == length(λs) "New interpolation coordinates θs::Vector, λs::Vector have to be of same length.
                                $n and $(length(λs)) provided."
    
    I = Interpolator(A,n)                   # generate Interpolator, containing geometry and work arrays
    update_locator!(θs,λs,I,unsafe=false)   # update location work arrays in I
    interpolate(A,I)
end

function interpolate(A::AbstractGrid,I::Interpolator{NF}) where NF
    @unpack npoints = I.locator             # number of points to interpolate onto
    As = zeros(NF,npoints)                  # preallocate: onto θs,λs interpolated values of A
    interpolate!(As,A,I)
end

function interpolate!(  As::Vector{NF},                 # Out: interpolated values
                        A::Grid,                        # gridded values to interpolate from
                        I::Interpolator{NF,Grid},       # geometry info and work arrays       
                        ) where {NF<:AbstractFloat,Grid<:AbstractGrid}

    @unpack ij_as,ij_bs,ij_cs,ij_ds,Δabs,Δcds,Δys = I.locator
    @unpack npoints = I.geometry
    @boundscheck length(As) == length(ij_as) || throw(BoundsError)

    A_northpole, A_southpole = average_on_poles(NF,A,I.geometry.rings)

    @boundscheck extrema_in(ij_as,-1,npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_bs,-1,npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_cs,-1,npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_ds,-1,npoints) || throw(BoundsError)

    @inbounds for (k,(ij_a,ij_b,ij_c,ij_d,Δab,Δcd,Δy)) in enumerate(zip(ij_as,ij_bs,ij_cs,ij_ds,Δabs,Δcds,Δys))

        # index ij=0,-1 indicates north,south pole
        a,b = ij_a ==  0 ? (A_northpole,A_northpole) : (A[ij_a],A[ij_b])
        c,d = ij_c == -1 ? (A_southpole,A_southpole) : (A[ij_c],A[ij_d])

        # weighted anvil-shaped average of a,b,c,d points around i point to interpolate on
        As[k] = anvil_average(a,b,c,d,Δab,Δcd,Δy)
    end

    return As
end

function update_locator!(   I::Interpolator{NF,Grid},   # geometry info and work arrays
                            θs::Vector,
                            λs::Vector;
                            unsafe::Bool=false,         # true to disable safety checks
                            ) where {NF<:AbstractFloat,Grid<:AbstractGrid}

    # find latitude ring indices corresponding to interpolation points
    @unpack latd = I.geometry
    @unpack js,Δys = I.locator
    find_rings!(js,Δys,θs,latd;unsafe)

    # find grid incides ij for top, bottom and left, right grid points around (θ,λ)
    find_grid_indices!(I,λs)
end

function find_rings!(   js::Vector{<:Integer},
                        Δys::Vector,            # distance fractions to ring further south
                        θs::Vector,
                        latd::Vector{NF};
                        unsafe::Bool=false) where NF
    
    if ~unsafe
        θmin,θmax = extrema(θs)
        @assert θmin >= -90 "Latitudes θs are expected to be within [-90˚,90˚]; θ=$(θmin)˚ given."
        @assert θmax <= 90 "Latitudes θs are expected to be within [-90˚,90˚]; θ=$(θmax)˚ given."
        
        @assert isdecreasing(latd) "Latitudes latd are expected to be strictly decreasing."
        @assert latd[1] == 90 "Latitudes latd are expected to contain 90˚C, the north pole."
        @assert latd[end] <= -90 "Latitudes latd are expected to contain -90˚, the south pole."
    end

    find_rings_unsafe!(js,Δys,θs,latd)
end

function find_rings_unsafe!(js::Vector{<:Integer},  # Out: vector of ring indices
                            Δys::Vector,            # distance fractions to ring further south
                            θs::Vector,             # latitudes of points to interpolate onto
                            latd::Vector{NF},       # latitudes of rings (90˚ to -90˚, strictly decreasing)
                            ) where {NF<:AbstractFloat}

    @boundscheck length(js) == length(θs) || throw(BoundsError)
    @boundscheck length(js) == length(Δys) || throw(BoundsError)

    # find first search for every θ in θs in latd but reuse index j from previous θ
    # as a predictor to start the search. Search walk is therefore either d=1 (southward)
    # or d=-1 (northward)
    
    j = 1                                       # starting ring, 0-based as latd contains poles
    @inbounds for (iθ,θf) in enumerate(θs)      # one latitude θ after another
        θ = convert(NF,θf)                      # convert to latd's format
        d,c = θ <= latd[j] ? (1,<=) : (-1,>)    # search direction d, d=1:south, d=-1:north, comparison c
        while c(θ,latd[j])                      # check whether ring j has been crossed in search direction
            j += d                              # walk in direction d
        end           
        j -= max(0,d)                           # so that [j,j+1) contains the point
        js[iθ] = j-1                            # convert back to 1-based indexed rings
        Δys[iθ] = (latd[j]-θ) / (latd[j]-latd[j+1])
    end
end

# for testing only
function find_rings(θs::Vector,latd::Vector{NF}) where NF
    js = Vector{Int}(undef,length(θs))
    Δys = Vector{NF}(undef,length(θs))
    find_rings!(js,Δys,θs,latd)
    return js,Δys
end

function find_grid_indices!(I::Interpolator,
                            λs::Vector)

    @unpack js, ij_as, ij_bs, ij_cs, ij_ds = I.locator
    @unpack Δabs, Δcds = I.locator
    @unpack nlons, rings, lon_offsets, nlat = I.geometry

    @inbounds for (k,(λf,j)) in enumerate(zip(λs,js))

        λ = convert(eltype(lon_offsets),λf)

        # NORTHERN POINTS a,b
        if j == 0           # a,b are at the north pole
            ij_as[k] = 0     # use 0 as north pole flag
            ij_bs[k] = 0
        else
            # get in-ring index i for a, the next grid point to the left
            # and b the next grid point to the right, such that
            # λ ∈ [a,b); while in most cases i_a + 1 = i_b, across 0˚E this is not the case
            i_a, i_b, Δ = find_lon_indices(λ,lon_offsets[j],nlons[j])
            ij_as[k] = rings[j][i_a]    # index ij for a
            ij_bs[k] = rings[j][i_b]    # index ij for b
            Δabs[k] = Δ                 # distance fraction of λ between a,b
        end

        # SOUTHERN POINTS c,d
        if j == nlat        # c,d are at the south pole
            ij_cs[k] = -1    # use -1 as south pole flag
            ij_ds[k] = -1
        else
            # as above but for one ring further down
            i_c, i_d, Δ = find_lon_indices(λ,lon_offsets[j+1],nlons[j+1])
            ij_cs[k] = rings[j+1][i_c]  # index ij for c
            ij_ds[k] = rings[j+1][i_d]  # index ij for d
            Δcds[k] = Δ                 # distance fraction of λ between c,d
        end
    end
end

function find_lon_indices(  λ::NF,
                            λ₀::NF,
                            nlon::Int) where {NF<:AbstractFloat}

    Δλ = 360/nlon           # longitude spacing
    ix = (λ-λ₀)/Δλ          # grid index i but with fractional part
    i = floor(Int,ix)       # 0-based grid index to the left
    Δ = ix-i                # distance fraction from i to i+1

    # [a,b] so that λ ∈ [λa,λb), i.e. a is the next grid point to the left, b to the right
    i_a = mod(i,nlon) + 1           # convert to 1-based index
    i_b = mod(i+1,nlon) + 1         # use mod for periodicity
    return i_a,i_b,Δ
end
    
average_on_poles(A::AbstractGrid{NF},rings::Vector{<:UnitRange}) where NF = average_on_poles(NF,A,rings)

function average_on_poles(  ::Type{NF},
                            A::AbstractGrid,
                            rings::Vector{<:UnitRange{<:Integer}}
                            ) where {NF<:AbstractFloat}
    
    A_northpole = mean(view(A,rings[1]))
    A_southpole = mean(view(A,rings[end]))
    return convert(NF,A_northpole), convert(NF,A_southpole)
end

"""
    r = anvil_average(a,b,c,d,Δab,Δcd,Δy)

The bilinear average of a,b,c,d which are values at grid points
in an anvil-shaped configuration at location x, which is denoted
by Δab,Δcd,Δy, the fraction of distances between a-b,c-d, and ab-cd,
respectively. Note that a,c and b,d do not necessarily share the same
longitude/x-coordinate. See schematic:

        0..............1    # fraction of distance Δab between a,b
        |<  Δab   >|

0^      a -------- o - b    # anvil-shaped average of a,b,c,d at location x
.Δy                |
.                  |
.v                 x 
.                  |
1         c ------ o ---- d

          |<  Δcd >|
          0...............1 # fraction of distance Δcd between c,d


^ fraction of distance Δy between a-b and c-d.
"""
function anvil_average( a::NF,      # top left value
                        b::NF,      # top right value
                        c::NF,      # bottom left value
                        d::NF,      # bottom right value
                        Δab::NF2,   # fraction of distance between a,b ∈ [0,1)
                        Δcd::NF2,   # fraction of distance between c,d ∈ [0,1)
                        Δy::NF2,    # fraction of distance between ab,cd ∈ [0,1)
                        ) where {NF,NF2}

    ab_average = a + (b-a)*convert(NF,Δab)      # a for Δab=1, b for Δab=1
    cd_average = c + (d-c)*convert(NF,Δcd)      # c for Δab=1, b for Δab=1
    abcd_average = ab_average + (cd_average-ab_average)*convert(NF,Δy)
    return abcd_average
end