abstract type AbstractGridGeometry end

"""Contains general precomputed arrays describing the grid.
$(TYPEDFIELDS)"""
struct GridGeometry{
    Grid,
    VectorType,
    VectorIntType,
} <: AbstractGridGeometry
    grid::Grid                  # grid, e.g. FullGaussianGrid

    nlat_half::Int              # number of latitude rings on one hemisphere (Eq. incl)
    nlat::Int                   # total number of latitude rings
    npoints::Int          # total number of grid points
    londs::VectorType           # longitudes of every grid point 0˚ to 360˚E
    latd::VectorType            # latitude of each ring, incl north pole 90˚N, ..., south pole -90˚N

    nlons::VectorIntType        # number of longitudinal points per ring
    lon_offsets::VectorType     # longitude offsets of first grid point per ring

    # rings, CPU copy of grid.rings
    rings::Vector{UnitRange{Int}} 
end

GridGeometry(field::AbstractField; kwargs...) = GridGeometry(field.grid; kwargs...)

"""
$(TYPEDSIGNATURES)          
Precomputed arrays describing the geometry of the Grid with resolution nlat_half.
Contains longitudes, latitudes of grid points, their ring index j and their
unravelled indices ij."""
function GridGeometry(
    grid::AbstractGrid;                                     # which grid to calculate the geometry for
    NF::Type{<:AbstractFloat} = DEFAULT_NF,                 # number format for the coordinates
)
    (; architecture, nlat_half) = grid
    nlat = get_nlat(grid)                       # total number of latitude rings
    npoints = get_npoints(grid)                 # total number of grid points

    # LATITUDES
    latd = get_latd(grid)                       # latitudes in degrees 90˚N ... -90˚N
    latd_poles = cat(90, latd, -90, dims=1)     # latd, but poles incl

    # Hack: use -90.00...1˚N instead of exactly -90˚N for the <=, > comparison
    # in find_rings! that way the last ring to the south pole can be an open
    # interval [a, b) with a the latitude of the last ring and b=90˚N as for
    # all other intervals between rings
    latd_poles[end] = latd_poles[end] - eps(latd_poles[end])

    # COORDINATES for every grid point in ring order
    londs, _ = get_londlatds(grid)              # in degrees [0˚...360˚E]                         

    # RINGS and LONGITUDE OFFSETS
    nlons = get_nlons(grid)                                 # number of longitude per ring, pole to pole
    lon_offsets = [londs[ring[1]] for ring in eachring(grid)]   # offset of the first point from 0˚E

    # vector type
    VectorType = array_type(architecture, NF, 1)
    VectorIntType = array_type(architecture, Int, 1)

    return GridGeometry{typeof(grid), VectorType, VectorIntType}(
        grid, nlat_half, nlat, npoints, londs, latd_poles, nlons, lon_offsets, Vector(grid.rings))
end

Base.show(io::IO,G::GridGeometry) = print(io,"GridGeometry for $(G.grid)")

"""Supertype of every Locator, which locates the indices on a grid to be used to perform an
interpolation. E.g. AnvilLocator uses a 4-point stencil for every new coordinate to interpolate
onto. Higher order stencils can be implemented by defining OtherLocator <: AbstractLocactor."""
abstract type AbstractLocator end

"""Contains arrays that locates grid points of a given field to be uses in an interpolation
and their weights. This Locator is a 4-point average in an anvil-shaped grid-point arrangement
between two latitude rings."""
@kwdef struct AnvilLocator{
    NF,
    VectorType,
    VectorIntType,
} <: AbstractLocator

    npoints_output::Int            # number of points to interpolate onto (length of following vectors)

    # to the coordinates respective indices
    js::VectorIntType       = zeros(Int, npoints_output)   # ring indices j such that [j, j+1) contains the point
    ij_as::VectorIntType    = zeros(Int, npoints_output)   # pixel index ij for top left point a on ring j
    ij_bs::VectorIntType    = zeros(Int, npoints_output)   # pixel index ij for top right point b on ring j
    ij_cs::VectorIntType    = zeros(Int, npoints_output)   # pixel index ij for bottom left point c on ring j+1
    ij_ds::VectorIntType    = zeros(Int, npoints_output)   # pixel index ij for bottom right point d on ring j+1

    # distances to adjacent grid points (i.e. the averaging weights)
    Δys::VectorType         = zeros(NF, npoints_output)    # distance fractions between rings
    Δabs::VectorType        = zeros(NF, npoints_output)    # distance fractions between a, b
    Δcds::VectorType        = zeros(NF, npoints_output)    # distance fractions between c, d
end

"""
$(TYPEDSIGNATURES)
Zero generator function for the 4-point average AnvilLocator. Use update_locator! to
update the grid indices used for interpolation and their weights. The number format
NF is the format used for the calculations within the interpolation, the input data
and/or output data formats may differ."""
function (::Type{L})(
    NF::Type{<:AbstractFloat},                                 # number format
    npoints::Integer; 
    architecture::AbstractArchitecture = DEFAULT_ARCHITECTURE(), # architecture to use
) where {L<:AbstractLocator}

    VectorType = array_type(architecture, NF, 1)
    VectorIntType = array_type(architecture, Int, 1)

    return L{NF, VectorType, VectorIntType}(;npoints_output=npoints)
end

# use Float32 as default for weights
(::Type{L})(npoints::Integer; kwargs...) where {L<:AbstractLocator} = L(DEFAULT_NF, npoints; kwargs...)

function Base.show(io::IO,L::AnvilLocator)
    println(io,"$(typeof(L))")
    print(io,"└ npoints_output::Int = $(L.npoints_output)")
end


"""
    abstract type AbstractInterpolator end

Supertype for Interpolators. Every Interpolator <: AbstractInterpolator is
expected to have two fields,
 - geometry, which describes the grid G to interpolate from
 - locator, which locates the indices on G and their weights to interpolate
    onto a new grid.
    
NF is the number format used to calculate the interpolation, which can be
different from the input data and/or the interpolated data on the new grid."""
abstract type AbstractInterpolator end

struct AnvilInterpolator{NF, Geometry, Locator} <: AbstractInterpolator
    geometry::Geometry
    locator::Locator
end

const DEFAULT_INTERPOLATOR = AnvilInterpolator
Base.eltype(::AnvilInterpolator{NF}) where NF = NF
grid_type(I::AnvilInterpolator) = typeof(I.geometry.grid)

function Base.show(io::IO, L::AnvilInterpolator{NF}) where NF
    println(io,"AnvilInterpolator{$NF} for $(L.geometry.grid)")
    print(io,"└ onto: $(L.locator.npoints_output) points")
end

# define to a <:AbstractInterpolator the corresponding Locator
Locator(::Type{<:AnvilInterpolator}) = AnvilLocator

# generator with NF default based on geometry and locator
function AnvilInterpolator(
    geometry::AbstractGridGeometry,
    locator::AbstractLocator;
    NF = DEFAULT_NF)
    return AnvilInterpolator{NF, typeof(geometry), typeof(locator)}(geometry, locator)
end

# generator from grid and npoints
function AnvilInterpolator(
    grid::AbstractGrid,
    npoints::Integer;       # number of points to interpolate onto
    NF::Type{<:AbstractFloat} = DEFAULT_NF,
)
    geometry = GridGeometry(grid)        # general coordinates and indices for grid
    locator = AnvilLocator(NF, npoints; architecture = grid.architecture)  # preallocate work arrays for interpolation

    # assemble geometry and locator to interpolator
    return AnvilInterpolator(geometry, locator; NF)
end

# generator that allocates a grid from Grid and nlat_half, TODO needed?
function (::Type{I})(   ::Type{Grid},
                        nlat_half::Integer,
                        npoints::Integer;
                        kwargs...
                        ) where {I<:AbstractInterpolator, Grid<:AbstractGrid}
    return I(Grid(nlat_half), npoints; kwargs...)
end

# general interpolator function that takes the Interpolator type as kwarg
function interpolator(
    grid::AbstractGrid,
    npoints::Integer;
    Interpolator::Type{<:AbstractInterpolator}=DEFAULT_INTERPOLATOR,
    kwargs...
)
    # note this leaves the interpolator.locator uninitialized
    return Interpolator(grid, npoints; kwargs...)
end

# generator based on two grids
function interpolator(
    grid_out::AbstractGrid,
    grid_in::AbstractGrid;
    kwargs...
)
    I = interpolator(grid_in, get_npoints(grid_out); kwargs...)
    londs, latds = get_londlatds(grid_out)
    londs = on_architecture(architecture(grid_out), londs)
    latds = on_architecture(architecture(grid_out), latds)

    update_locator!(I, londs, latds, unsafe=false)
    return I
end

# generator based on two fields
interpolator(Aout::AbstractField, A::AbstractField; kwargs...) = 
    interpolator(Aout.grid, A.grid; kwargs...)
    
## FUNCTIONS
# interpolate into one point, wrap them into vector
function interpolate(lond::Real, latd::Real, A::AbstractField2D; kwargs...)
    Ai = interpolate([lond], [latd], A; kwargs...)
    return Ai[1]
end

# interpolate while creating an interpolator on the fly
function interpolate(
    londs::AbstractVector{NF},  # longitudes to interpolate into (0˚...360˚E)
    latds::AbstractVector{NF},  # latitudes to interpolate onto (90˚N...-90˚N)
    A::AbstractField2D;         # gridded field to interpolate from
    kwargs...
) where NF                      # number format used for interpolation
    npoints = length(latds)
    npoints == length(londs) || DimensionMismatch("New interpolation coordinates londs, latds have to be of same length.
                                $n and $(length(londs)) provided.")
    
    I = interpolator(A.grid, npoints; kwargs...)    # generate Interpolator, containing geometry and work arrays
    update_locator!(I, londs, latds, unsafe=false)  # update location work arrays in I
    interpolate(A, I)
end

# interpolate using an existing interpolator (needs to be initialized with update_locator!)
function interpolate(
    A::AbstractField2D{NF},     # field to interpolate 
    I::AbstractInterpolator;    # indices in I are assumed to be calculated already!
    kwargs...
) where NF
    (; npoints ) = I.locator                                         # number of points to interpolate onto
    Aout = array_type(architecture(A), NF, 1)(undef, npoints)    # preallocate: onto θs, λs interpolated values of A
    _interpolate!(Aout, A.data, I, architecture=architecture(A)) # perform interpolation, store in Aout
end

# the actual interpolation function
function _interpolate!(
    Aout,               # Out: interpolated values
    A,                  # gridded values to interpolate from
    interpolator::AnvilInterpolator,    # geometry info and work arrays 
    architecture::AbstractArchitecture      
)
    (; npoints_output, ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys) = interpolator.locator
    (; npoints, rings) = interpolator.geometry
    
    # 1) Aout's length must match the interpolator
    # 2) input A must match the interpolator's geometry points (do not check grids for view support)
    @boundscheck length(Aout) == length(ij_as) || throw(DimensionMismatchArray(Aout, ij_as))
    @boundscheck length(A) == npoints ||
        throw(DimensionMismatch("Interpolator ($npoints points) mismatches input grid ($(length(A)) points)."))

    A_northpole, A_southpole = average_on_poles(A, rings)

    #TODO ij_cs, ij_ds shouldn't be 0...
    @boundscheck extrema_in(ij_as,  0, npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_bs,  0, npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_cs, -1, npoints) || throw(BoundsError)
    @boundscheck extrema_in(ij_ds, -1, npoints) || throw(BoundsError)

    launch!(architecture,
    :linear,
    (npoints_output,),
    _interpolate_kernel!,
        Aout,
        A,
        ij_as,
        ij_bs,
        ij_cs,
        ij_ds,
        Δabs,
        Δcds,
        Δys,
        A_northpole,
        A_southpole,
    )
    synchronize(architecture)

    return Aout
end

@kernel inbounds=true function _interpolate_kernel!(
    Aout,               # Out: interpolated values
    A,                  # gridded values to interpolate from
    @Const(ij_as),      # indices of A to interpolate from
    @Const(ij_bs),      # indices of A to interpolate from
    @Const(ij_cs),      # indices of A to interpolate from
    @Const(ij_ds),      # indices of A to interpolate from
    @Const(Δabs),       # weights of A to interpolate from
    @Const(Δcds),       # weights of A to interpolate from
    @Const(Δys),        # weights of A to interpolate from
    @Const(A_northpole),
    @Const(A_southpole),
)
    
    k = @index(Global, Linear)

    a, b = ij_as[k] ==  0 ? (A_northpole, A_northpole) : (A[ij_as[k]], A[ij_bs[k]])
    c, d = ij_cs[k] == -1 ? (A_southpole, A_southpole) : (A[ij_cs[k]], A[ij_ds[k]])
    
    Aout[k] = anvil_average(a, b, c, d, Δabs[k], Δcds[k], Δys[k])
end

# version for 3D+ fields GPU 
function interpolate!(
    Aout::AbstractField,
    A::AbstractField,
    interpolator::AbstractInterpolator,
) 
    fields_match(Aout, A) && return copyto!(Aout.data, A.data)
    @assert ismatching(architecture(A), A_out) "Interpolation is only supported between fields on the same architecture."
    _interpolate!(Aout.data, A.data, interpolator, architecture(A))
end

# interpolate while creating an interpolator on the fly
function interpolate!(
    Aout::AbstractField,
    A::AbstractField;
    kwargs...
)
    # if fields match just copy data over (eltypes might differ)
    fields_match(Aout, A) && return copyto!(Aout.data, A.data)
    I = interpolator(Aout, A; kwargs...)    # create interpolator instance from field A to Aout
    interpolate!(Aout, A, I)                # perform interpolation returns 
    return Aout                             # return the field wrapped around the interpolated data
end

# create grid on the fly
interpolate(Grid::Type{<:AbstractGrid}, nlat_half::Integer, A::AbstractField; kwargs...) =
    interpolate(Grid(nlat_half), A; kwargs...)

# create field from grid on the fly
function interpolate(
    grid::AbstractGrid,
    A::AbstractField;
    kwargs...
)
    I = interpolator(grid, A.grid; kwargs...)
    Aout = Field(grid, size(A)[2:end]...)
    interpolate!(Aout, A, I)    # returns an array
    return Aout                 # returns the field wrapped around that array
end

# if only the grid type is provided, create a grid with nlat_half and architecture from the input field
interpolate(Grid::Type{<:AbstractGrid}, A::AbstractField; kwargs...) = interpolate(Grid(A.grid.nlat_half, architecture(A)), A; kwargs...)

function update_locator!(
    I::AbstractInterpolator,    # GridGeometry and Locator
    λs::AbstractVector,         # longitudes to interpolate onto
    θs::AbstractVector;         # latitudes to interpolate onto
    unsafe::Bool=false,         # true to disable safety checks
)
    # find latitude ring indices corresponding to interpolation points
    (; latd ) = I.geometry                  # latitudes of rings including north and south pole
    (; js, Δys ) = I.locator                # to be updated: ring indices js, and meridional weights Δys
    find_rings!(js, Δys, θs, latd; unsafe, architecture=I.geometry.grid.architecture)  # next ring at or north of θ

    # find grid incides ij for top, bottom and left, right grid points around (θ, λ)
    find_grid_indices!(I, λs, I.geometry.grid.architecture)               # next points left and right of λ on rings north and south
end

function update_locator!(I::AbstractInterpolator, A::AbstractField; kwargs...)
    londs, latds = get_londlatds(A.grid)
    londs = on_architecture(architecture(A), londs)
    latds = on_architecture(architecture(A), latds)

    update_locator!(I, londs, latds; kwargs...)
end

function find_rings!(   js::AbstractVector{<:Integer},  # Out: ring indices j
                        Δys::AbstractVector,            # Out: distance fractions to ring further south
                        θs::AbstractVector,             # latitudes to interpolate onto
                        latd::AbstractVector;           # latitudes of the rings on the original grid
                        unsafe::Bool=false,             # skip safety checks when true
                        architecture::AbstractArchitecture=architecture(js))     
    
    if ~unsafe
        θmin, θmax = extrema(θs)
        @assert θmin >= -90 "Latitudes θs are expected to be within [-90˚, 90˚]; θ=$(θmin)˚ given."
        @assert θmax <= 90 "Latitudes θs are expected to be within [-90˚, 90˚]; θ=$(θmax)˚ given."
        
        #TODO: currently we only check the latitudes of the grid we interpolate onto 
        #TODO: as we only allow instances of Field to be the original grid, the latitudes below 
        #TODO: should be okay anyway. Checking it on GPU is a bit more annoying...

        #@assert isdecreasing(latd) "Latitudes latd are expected to be strictly decreasing."
        #@assert latd[1] == 90 "Latitudes latd are expected to contain 90˚N, the north pole."

        # Hack: for intervals between rings to be one-sided open [j, j+1) the last element in
        # latd has to be prevfloat(-90) for the <=, > comparisons
        #ϵ = eps(latd[end])  
        #@assert latd[end] == -90-ϵ "Latitudes latd are expected to contain -90˚, the south pole."
    end

    find_rings_unsafe!(js, Δys, θs, latd, architecture)
end

DimensionMismatchArray(a::AbstractArray, b::AbstractArray) =
    DimensionMismatch("Arrays have different dimensions: $(size(a)) vs $(size(b))")
DimensionMismatchArray(a::AbstractArray, bs::AbstractArray...) =
    DimensionMismatch("Arrays have different dimensions: $(size(a)) vs $(map(size, bs))")

# CPU version of find_rings_unsafe!
function find_rings_unsafe!(js::AbstractVector{<:Integer},  # Out: vector of ring indices
                            Δys::AbstractVector,            # distance fractions to ring further south
                            θs::AbstractVector,             # latitudes of points to interpolate onto
                            latd::AbstractVector{NF},       # latitudes of rings (90˚ to -90˚, strictly decreasing)
                            architecture::CPU) where {NF<:AbstractFloat}

    @boundscheck length(js) == length(θs) || throw(DimensionMismatchArray(js, θs))
    @boundscheck length(js) == length(Δys) || throw(DimensionMismatchArray(js, Δys))

    # find first search for every θ in θs in latd but reuse index j from previous θ
    # as a predictor to start the search. Search walk is therefore either d=1 (southward)
    # or d=-1 (northward)
    
    j = 1                                       # starting ring, 0-based as latd contains poles
    @inbounds for (iθ, θf) in enumerate(θs)     # one latitude θ after another
        θ = convert(NF, θf)                     # convert to latd's format
        d, c = θ <= latd[j] ? (1, <=) : (-1, >) # search direction d, d=1:south, d=-1:north, comparison c
        while c(θ, latd[j])                     # check whether ring j has been crossed in search direction
            j += d                              # walk in direction d
        end           
        j -= max(0, d)                          # so that [j, j+1) contains the point
        js[iθ] = j-1                            # convert back to 1-based indexed rings
        Δys[iθ] = (latd[j]-θ) / (latd[j]-latd[j+1])
    end
end

# GPU kernel for finding ring indices
@kernel inbounds=true function find_rings_kernel!(
    js,                # Out: ring indices j
    Δys,               # Out: distance fractions to ring further south
    @Const(θs),        # latitudes to interpolate onto
    @Const(latd)       # latitudes of the rings on the original grid
)
    k = @index(Global, Linear)
    
    # Binary search to find the ring index for this latitude
    θ = θs[k]
    
    # Initialize search bounds
    j_low = 1
    j_high = length(latd) - 1
    
    # Binary search for the ring that contains this latitude
    j = 1  # Default starting point
    
    # Binary search is more efficient on GPU than linear search
    while j_low <= j_high
        j_mid = (j_low + j_high) ÷ 2
        
        if θ <= latd[j_mid] && (j_mid == length(latd) - 1 || θ > latd[j_mid + 1])
            # Found the correct ring
            j = j_mid
            break
        elseif θ > latd[j_mid]
            j_high = j_mid - 1
        else
            j_low = j_mid + 1
        end
    end
    
    # Store results
    js[k] = j - 1  # Convert to 1-based indexed rings
    Δys[k] = (latd[j] - θ) / (latd[j] - latd[j+1])
end

# GPU version of find_rings_unsafe! using KernelAbstractions
function find_rings_unsafe!(js::AbstractArray{<:Integer},  # Out: vector of ring indices
                           Δys::AbstractArray,            # distance fractions to ring further south
                           θs::AbstractArray,             # latitudes of points to interpolate onto
                           latd::AbstractArray{NF},       # latitudes of rings (90˚ to -90˚, strictly decreasing)
                           architecture::GPU) where {NF<:AbstractFloat}
    
    @boundscheck length(js) == length(θs) || throw(DimensionMismatchArray(js, θs))
    @boundscheck length(js) == length(Δys) || throw(DimensionMismatchArray(js, Δys))
        
    launch!(
        architecture,
        :linear,
        size(js),
        find_rings_kernel!,
        js,
        Δys,
        θs,
        latd
    )
    synchronize(architecture)
end

# for testing only
function find_rings(θs::AbstractVector, latd::AbstractVector{NF}) where NF
    js = Vector{Int}(undef, length(θs))
    Δys = Vector{NF}(undef, length(θs))
    find_rings!(js, Δys, θs, latd)
    return js, Δys
end

# CPU version of find_grid_indices!
function find_grid_indices!(I::AnvilInterpolator,       # update indices arrays
                            λs::AbstractVector, 
                            architecture::CPU)         # based on new longitudes λ

    (; js, ij_as, ij_bs, ij_cs, ij_ds ) = I.locator
    (; Δabs, Δcds ) = I.locator
    (; nlons, lon_offsets, nlat ) = I.geometry
    (; rings ) = I.geometry.grid

    @inbounds for (k, (λf, j)) in enumerate(zip(λs, js))

        λ = convert(eltype(lon_offsets), λf)

        # NORTHERN POINTS a, b
        if j == 0               # a, b are at the north pole
            ij_as[k] = 0        # use 0 as north pole flag
            ij_bs[k] = 0
        else
            # get in-ring index i for a, the next grid point to the left
            # and b the next grid point to the right, such that
            # λ ∈ [a, b); while in most cases i_a + 1 = i_b, across 0˚E this is not the case
            i_a, i_b, Δ = find_lon_indices(λ, lon_offsets[j], nlons[j])
            ij_as[k] = rings[j][i_a]    # index ij for a
            ij_bs[k] = rings[j][i_b]    # index ij for b
            Δabs[k] = Δ                 # distance fraction of λ between a, b
        end

        # SOUTHERN POINTS c, d
        if j == nlat                    # c, d are at the south pole
            ij_cs[k] = -1               # use -1 as south pole flag
            ij_ds[k] = -1
        else
            # as above but for one ring further down
            i_c, i_d, Δ = find_lon_indices(λ, lon_offsets[j+1], nlons[j+1])
            ij_cs[k] = rings[j+1][i_c]  # index ij for c
            ij_ds[k] = rings[j+1][i_d]  # index ij for d
            Δcds[k] = Δ                 # distance fraction of λ between c, d
        end
    end
end

# GPU kernel for finding grid indices
@kernel inbounds=true function find_grid_indices_kernel!(
    @Const(js),            # ring indices j
    ij_as, ij_bs,          # northern point indices
    ij_cs, ij_ds,          # southern point indices
    Δabs, Δcds,            # distance fractions
    @Const(λs),            # longitudes to interpolate onto
    @Const(lon_offsets),   # longitude offsets for each ring
    @Const(nlons),         # number of longitude points per ring
    @Const(nlat),          # number of latitude rings
    @Const(rings)          # ring indices
)
    k = @index(Global, Linear)
    
    j = js[k]
    λ = λs[k]
    
    # NORTHERN POINTS a, b
    if j == 0               # a, b are at the north pole
        ij_as[k] = 0        # use 0 as north pole flag
        ij_bs[k] = 0
    else
        # get in-ring index i for a, the next grid point to the left
        # and b the next grid point to the right, such that
        # λ ∈ [a, b); while in most cases i_a + 1 = i_b, across 0˚E this is not the case
        i_a, i_b, Δ = find_lon_indices(λ, lon_offsets[j], nlons[j])
        ij_as[k] = rings[j][i_a]    # index ij for a
        ij_bs[k] = rings[j][i_b]    # index ij for b
        Δabs[k] = Δ                 # distance fraction of λ between a, b
    end

    # SOUTHERN POINTS c, d
    if j == nlat                    # c, d are at the south pole
        ij_cs[k] = -1               # use -1 as south pole flag
        ij_ds[k] = -1
    else
        # as above but for one ring further down
        i_c, i_d, Δ = find_lon_indices(λ, lon_offsets[j+1], nlons[j+1])
        ij_cs[k] = rings[j+1][i_c]  # index ij for c
        ij_ds[k] = rings[j+1][i_d]  # index ij for d
        Δcds[k] = Δ                 # distance fraction of λ between c, d
    end
end

# GPU version of find_grid_indices! using KernelAbstractions
function find_grid_indices!(I::AnvilInterpolator,       # update indices arrays
                           λs::AbstractArray,           # based on new longitudes λ
                           architecture::GPU)
                           
    (; js, ij_as, ij_bs, ij_cs, ij_ds ) = I.locator
    (; Δabs, Δcds ) = I.locator
    (; nlons, lon_offsets, nlat ) = I.geometry
    (; rings ) = I.geometry.grid
    
    # Convert λs to the same type as lon_offsets if needed
    λs_converted = convert.(eltype(lon_offsets), λs)
    
    launch!(
        architecture,
        :linear,
        size(js),
        find_grid_indices_kernel!,
        js,
        ij_as, ij_bs,
        ij_cs, ij_ds,
        Δabs, Δcds,
        λs_converted,
        lon_offsets,
        nlons,
        nlat,
        rings
    )
    synchronize(architecture)
end

@inline function find_lon_indices(λ::NF,     # longitude to find incides for (0˚...360˚E)
                            λ₀::NF,     # offset of the first longitude point on ring
                            nlon::Int   # number of longitude points on ring
                            ) where {NF<:AbstractFloat}

    Δλ = convert(NF, 360)/nlon          # longitude spacing
    ix = (λ-λ₀)/Δλ                      # grid index i but with fractional part
    i = floor(Int, ix)                  # 0-based grid index to the left
    Δ = ix-i                            # distance fraction from i to i+1

    # λ ∈ [λa, λb), i.e. a is the next grid point to the left, b to the right
    i_a = mod(i, nlon) + 1              # convert to 1-based index
    i_b = mod(i+1, nlon) + 1            # use mod for periodicity
    return i_a, i_b, Δ
end

"""
$(TYPEDSIGNATURES)
Computes the average at the North and South pole from a given grid `A` and it's precomputed
ring indices `rings`. The North pole average is an equally weighted average of all grid points
on the northern-most ring. Similar for the South pole."""
@inline function average_on_poles(A::AbstractVector{NF}, rings) where {NF<:AbstractFloat} 
    # TODO: doing the computation below with views causes scalarindexing on GPUs
    A_northpole = mean(A[rings[1]])     # average of all grid points around the north pole
    A_southpole = mean(A[rings[end]])   # same for south pole
    return A_northpole, A_southpole
end
    
"""
$(TYPEDSIGNATURES)
Method for `A::Abstract{T<:Integer}` which rounds the averaged values
to return the same number format `NF`."""
@inline function average_on_poles(A::AbstractVector{NF}, rings) where {NF<:Integer}
    # TODO: doing the computation below with views causes scalarindexing on GPUs
    A_northpole = mean(A[rings[1]])    # average of all grid points around the north pole
    A_southpole = mean(A[rings[end]])  # same for south pole
    return round(NF, A_northpole), round(NF, A_southpole)
end

"""
$(TYPEDSIGNATURES)
The bilinear average of a, b, c, d which are values at grid points
in an anvil-shaped configuration at location x, which is denoted
by Δab, Δcd, Δy, the fraction of distances between a-b, c-d, and ab-cd,
respectively. Note that a, c and b, d do not necessarily share the same
longitude/x-coordinate. See schematic:
```
            0..............1    # fraction of distance Δab between a, b
            |<  Δab   >|

    0^      a -------- o - b    # anvil-shaped average of a, b, c, d at location x
    .Δy                |
    .                  |
    .v                 x 
    .                  |
    1         c ------ o ---- d

              |<  Δcd >|
              0...............1 # fraction of distance Δcd between c, d
```
^ fraction of distance Δy between a-b and c-d."""
@inline function anvil_average(
    a,      # top left value
    b,      # top right value
    c,      # bottom left value
    d,      # bottom right value
    Δab,    # fraction of distance between a, b ∈ [0, 1)
    Δcd,    # fraction of distance between c, d ∈ [0, 1)
    Δy,     # fraction of distance between ab, cd ∈ [0, 1)
    )

    # the type of the weights is ::Real, but the following is written such that
    # always NF (the type of the data values a, b, c, d) is returned
    ab_average = a + (b-a)*Δab      # a for Δab=0, b for Δab=1
    cd_average = c + (d-c)*Δcd      # c for Δab=0, b for Δab=1
    abcd_average = ab_average + (cd_average-ab_average)*Δy
    return abcd_average
end

"""
$TYPEDSIGNATURES
Averages all grid points in `input` that are within one grid cell of
`output` with coslat-weighting. The output grid cell boundaries 
are assumed to be rectangles spanning half way to adjacent longitude
and latitude points."""
function grid_cell_average!(
    output::AbstractField,
    input::AbstractFullField)
        
    # for i, j indexing
    input_matrix = reshape(input.data, :, get_nlat(input))

    # input grid coordinates, all in radians 0:π for colat, 0:2π for lon
    colat_in = get_colat(input)
    coslat = sin.(colat_in)    # cos(lat) = sin(colat)
    lon_in = get_lon(input)
    nlon_in = length(lon_in)

    # output grid coordinates, append -π, 2π to have grid points
    # towards the poles definitely included
    colat_out = vcat(-π, get_colat(output), 2π)
    lons_out, _ = get_lonlats(output)

    rings = eachring(output)

    for (j, ring) in enumerate(rings)
        Δϕ = 2π/length(ring)        # longitude spacing on this ring

        # indices for lat_out are shifted as north and south pole are included
        θ0 = (colat_out[j]   + colat_out[j+1])/2    # northern edge
        θ1 = (colat_out[j+1] + colat_out[j+2])/2    # southern edge

        # matrix indices for input grid that lie in output grid cell
        j0 = findfirst(θ -> θ >= θ0, colat_in)
        j1 = findlast( θ -> θ <  θ1, colat_in)

        for ij in ring
            ϕ0 = mod(lons_out[ij] - Δϕ/2, 2π)       # western edge
            ϕ1 = ϕ0 + Δϕ                            # eastern edge
            # the last line does not require a mod and in fact throws
            # an error if, as long as the offset from prime meridian
            # is at most Δϕ/2 (true for HEALPixGrids, offset 0 for
            # Gaussian and Clenshaw grids)

            # matrix indices for input grid that lie in output grid cell
            a = findlast( ϕ -> ϕ <  ϕ0, lon_in)     # western edge
            b = findfirst(ϕ -> ϕ >= ϕ1, lon_in)     # eastern edge

            # map around prime meridian if coordinates outside of range
            a = isnothing(a) ? nlon_in : a
            b = isnothing(b) ? 1 : b
            
            # in most cases we will have 1 <= a < b <=n, then loop i in a:b (b:a isn't looping)
            # however at the edges we have a < 1 or n < b the mod turns this into
            # 1 <= b < a <= n, for which we want to loop 1:b AND a:n
            ab, ba = b < a ? (1:b, a:nlon_in) : (a:b, b:a)

            sum_of_weights = 0
            for j′ in j0:j1
                    
                for i in ab
                    w = coslat[j′]
                    sum_of_weights += w
                    output[ij] += w*input_matrix[i, j′]
                end
                
                for i in ba
                    w = coslat[j′]
                    sum_of_weights += w
                    output[ij] += w*input_matrix[i, j′]
                end
            end
            output[ij] /= sum_of_weights
        end
    end
    
    return output
end


"""
$TYPEDSIGNATURES
Averages all grid points in `input` that are within one grid cell of
`output` with coslat-weighting. The output grid cell boundaries 
are assumed to be rectangles spanning half way to adjacent longitude
and latitude points."""
function grid_cell_average(Grid::Type{<:AbstractGrid}, nlat_half::Integer, input::AbstractFullGrid)
    grid = Grid(nlat_half)
    output = Field(grid)
    grid_cell_average!(output, input)
    return output
end