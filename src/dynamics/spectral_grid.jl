const DEFAULT_NF = Float32
const DEFAULT_MODEL = PrimitiveDry
const DEFAULT_GRID = OctahedralGaussianGrid

"""
Defines the horizontal spectral resolution and corresponding grid and the
vertical coordinate for SpeedyWeather.jl. Options are
$(TYPEDFIELDS)

`nlat_half` and `npoints` should not be chosen but are derived from `trunc`,
`Grid` and `dealiasing`."""
Base.@kwdef struct SpectralGrid
    "number format used throughout the model"
    NF::Type{<:AbstractFloat} = DEFAULT_NF

    # HORIZONTAL
    "horizontal resolution as the maximum degree of spherical harmonics"
    trunc::Int = 31

    "horizontal grid used for calculations in grid-point space"
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID

    "how to match spectral with grid resolution: dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64 = 2

    "radius of the sphere [m]"
    radius::Float64 = 6.371e6

    # SIZE OF GRID from trunc, Grid, dealiasing:
    "number of latitude rings on one hemisphere (Equator incl)"
    nlat_half::Int = SpeedyTransforms.get_nlat_half(trunc,dealiasing)

    "total number of grid points in the horizontal"
    npoints::Int = RingGrids.get_npoints(Grid,nlat_half)

    # VERTICAL
    "number of vertical levels"
    nlev::Int = 8

    "coordinates used to discretize the vertical"
    vertical_coordinates::VerticalCoordinates = SigmaCoordinates(;nlev)

    # make sure nlev and vertical_coordinates.nlev match
    function SpectralGrid(NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,nlev,vertical_coordinates)
        if nlev == vertical_coordinates.nlev
            return new(NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,
                    nlev,vertical_coordinates)
        else    # use nlev from vert_coords:
            return new(NF,trunc,Grid,dealiasing,radius,nlat_half,npoints,
                    vertical_coordinates.nlev,vertical_coordinates)
        end
    end
end

# generator functions
SpectralGrid(NF::Type{<:AbstractFloat};kwargs...) = SpectralGrid(;NF,kwargs...)
SpectralGrid(Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(;Grid,kwargs...)
SpectralGrid(NF::Type{<:AbstractFloat},Grid::Type{<:AbstractGrid};kwargs...) = SpectralGrid(;NF,Grid,kwargs...)

function Base.show(io::IO,SG::SpectralGrid)
    (;NF,trunc,Grid,radius,nlat_half,npoints,nlev,vertical_coordinates) = SG
    # truncation = if dealiasing < 2 "linear" elseif dealiasing < 3 "quadratic" else "cubic" end
    
    # resolution information
    res_ave = sqrt(4π*radius^2/npoints)/1000  # in [km]
    res_eq_x = 2π*radius/RingGrids.get_nlon_max(Grid,nlat_half)/1000
    lat = get_lat(Grid,nlat_half)
    res_eq_y = (lat[nlat_half] - lat[nlat_half+1])*radius/1000

    s(x) = @sprintf("%.3g",x)

    println(io,"$(typeof(SG)):")
    println(io,"├ Spectral:   T$trunc LowerTriangularMatrix{Complex{$NF}}, radius = $radius m")
    println(io,"├ Grid:       $(get_nlat(Grid,nlat_half))-ring $Grid{$NF}, $npoints grid points")
    println(io,"├ Resolution: $(s(res_ave))km (average), $(s(res_eq_x))km × $(s(res_eq_y))km (Equator)")
      print(io,"└ Vertical:   $nlev-level $(typeof(vertical_coordinates))")
end

"""
$(TYPEDSIGNATURES)
Construct Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Geometry{NF<:AbstractFloat} <: AbstractGeometry{NF}     # NF: Number format

    "SpectralGrid that defines spectral and grid resolution"
    spectral_grid::SpectralGrid

    "grid of the dynamical core"
    Grid::Type{<:AbstractGrid} = spectral_grid.Grid

    "resolution parameter nlat_half of Grid, # of latitudes on one hemisphere (incl Equator)"
    nlat_half::Int = spectral_grid.nlat_half      


    # GRID-POINT SPACE
    "maximum number of longitudes (at/around Equator)"
    nlon_max::Int = get_nlon_max(Grid,nlat_half)

    "=nlon_max, same (used for compatibility), TODO: still needed?"
    nlon::Int = nlon_max

    "number of latitude rings"
    nlat::Int = get_nlat(Grid,nlat_half)

    "number of vertical levels"
    nlev::Int = spectral_grid.nlev

    "total number of grid points"
    npoints::Int = spectral_grid.npoints

    "Planet's radius [m]"
    radius::NF = spectral_grid.radius    


    # ARRAYS OF LANGITUDES/LONGITUDES
    "array of colatitudes in radians (0...π)"
    colat::Vector{Float64} = get_colat(Grid,nlat_half)

    "array of latitudes in degrees (90˚...-90˚)"
    latd::Vector{Float64} = get_latd(Grid,nlat_half)

    "array of longitudes in degrees (0...360˚), empty for non-full grids"
    lond::Vector{Float64} = get_lond(Grid,nlat_half)

    "longitude (-180˚...180˚) for each grid point in ring order"
    londs::Vector{NF} = get_latdlonds(Grid,nlat_half)[2]
    
    "latitude (-90˚...˚90) for each grid point in ring order"
    latds::Vector{NF} = get_latdlonds(Grid,nlat_half)[1]

    "sin of latitudes"
    sinlat::Vector{NF} = sind.(latd)
    
    "cos of latitudes"
    coslat::Vector{NF} = cosd.(latd)
    
    "= 1/cos(lat)"
    coslat⁻¹::Vector{NF} = 1 ./ coslat

    "= cos²(lat)"
    coslat²::Vector{NF} = coslat.^2

    "# = 1/cos²(lat)"
    coslat⁻²::Vector{NF} = 1 ./ coslat²            


    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::Vector{NF} = spectral_grid.vertical_coordinates.σ_half

    "σ at full levels, σₖ"
    σ_levels_full::Vector{NF} = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])  
    
    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::Vector{NF} = σ_levels_half[2:end] - σ_levels_half[1:end-1]      

    "log of σ at full levels, include surface (σ=1) as last element"
    ln_σ_levels_full::Vector{NF} = log.(vcat(σ_levels_full,1))
end

"""
$(TYPEDSIGNATURES)
Generator function for `Geometry` struct based on `spectral_grid`."""
function Geometry(spectral_grid::SpectralGrid)
    return Geometry{spectral_grid.NF}(;spectral_grid)
end

function Base.show(io::IO,G::Geometry)
    print(io,"$(typeof(G)) for $(G.spectral_grid)")
end

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct pulling in parameters from a SpectralGrid struct."""
function SpeedyTransforms.SpectralTransform(spectral_grid::SpectralGrid;
                                            recompute_legendre::Bool = false,
                                            one_more_degree::Bool = true,
                                            kwargs...)
    (;NF, Grid, trunc, dealiasing) = spectral_grid
    return SpectralTransform(NF,Grid,trunc+one_more_degree,trunc;recompute_legendre,dealiasing,kwargs...)
end

