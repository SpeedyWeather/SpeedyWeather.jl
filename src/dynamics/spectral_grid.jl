abstract type AbstractSpectralGrid end
abstract type AbstractGeometry end

const DEFAULT_NF = Float32
const DEFAULT_MODEL = PrimitiveDry
const DEFAULT_GRID = OctahedralGaussianGrid
const DEFAULT_RADIUS = 6.371e6
const DEFAULT_TRUNC = 31
const DEFAULT_NLEV = 8

export SpectralGrid

"""
Defines the horizontal spectral resolution and corresponding grid and the
vertical coordinate for SpeedyWeather.jl. Options are
$(TYPEDFIELDS)

`nlat_half` and `npoints` should not be chosen but are derived from `trunc`,
`Grid` and `dealiasing`."""
Base.@kwdef struct SpectralGrid <: AbstractSpectralGrid
    "number format used throughout the model"
    NF::Type{<:AbstractFloat} = DEFAULT_NF

    # HORIZONTAL
    "horizontal resolution as the maximum degree of spherical harmonics"
    trunc::Int = DEFAULT_TRUNC

    "horizontal grid used for calculations in grid-point space"
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID

    "how to match spectral with grid resolution: dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64 = 2

    "radius of the sphere [m]"
    radius::Float64 = DEFAULT_RADIUS

    # SIZE OF GRID from trunc, Grid, dealiasing:
    "number of latitude rings on one hemisphere (Equator incl)"
    nlat_half::Int = SpeedyTransforms.get_nlat_half(trunc,dealiasing)

    "total number of grid points in the horizontal"
    npoints::Int = RingGrids.get_npoints(Grid,nlat_half)

    # VERTICAL
    "number of vertical levels"
    nlev::Int = DEFAULT_NLEV

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

    s(x) = x > 1000 ? @sprintf("%i",x) : @sprintf("%.3g",x)

    println(io,"$(typeof(SG)):")
    println(io,"├ Spectral:   T$trunc LowerTriangularMatrix{Complex{$NF}}, radius = $radius m")
    println(io,"├ Grid:       $(get_nlat(Grid,nlat_half))-ring $Grid{$NF}, $npoints grid points")
    println(io,"├ Resolution: $(s(res_ave))km (average), $(s(res_eq_x))km × $(s(res_eq_y))km (Equator)")
      print(io,"└ Vertical:   $nlev-level $(typeof(vertical_coordinates))")
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

