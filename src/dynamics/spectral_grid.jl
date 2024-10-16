abstract type AbstractSpectralGrid end

# computing
const DEFAULT_NF = Float32
const DEFAULT_DEVICE = CPU()
const DEFAULT_ARRAYTYPE = Array

# numerics
const DEFAULT_GRID = OctahedralGaussianGrid
const DEFAULT_RADIUS = 6.371e6
const DEFAULT_TRUNC = 31
const DEFAULT_NLAYERS = 8

export SpectralGrid

"""
Defines the horizontal spectral resolution and corresponding grid and the
vertical coordinate for SpeedyWeather.jl. Options are
$(TYPEDFIELDS)

`nlat_half` and `npoints` should not be chosen but are derived from `trunc`,
`Grid` and `dealiasing`."""
@kwdef struct SpectralGrid <: AbstractSpectralGrid
    "[OPTION] number format used throughout the model"
    NF::Type{<:AbstractFloat} = DEFAULT_NF

    "[OPTION] device archictecture to run on"
    device::AbstractDevice = DEFAULT_DEVICE

    "[OPTION] array type to use for all variables"
    ArrayType::Type{<:AbstractArray} = default_array_type(device)

    # HORIZONTAL
    "[OPTION] horizontal resolution as the maximum degree of spherical harmonics"
    trunc::Int = DEFAULT_TRUNC

    "[OPTION] horizontal grid used for calculations in grid-point space"
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID

    "[OPTION] how to match spectral with grid resolution: dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64 = 2

    "[OPTION] radius of the sphere [m]"
    radius::Float64 = DEFAULT_RADIUS

    "[OPTION] number of particles for particle advection [1]"
    nparticles::Int = 0

    # SIZE OF GRID from trunc, Grid, dealiasing:
    "number of latitude rings on one hemisphere (Equator incl)"
    nlat_half::Int = SpeedyTransforms.get_nlat_half(trunc, dealiasing)

    "number of latitude rings on both hemispheres"
    nlat::Int = RingGrids.get_nlat(Grid, nlat_half)

    "total number of grid points in the horizontal"
    npoints::Int = RingGrids.get_npoints(Grid, nlat_half)

    # VERTICAL
    "[OPTION] number of vertical levels"
    nlayers::Int = DEFAULT_NLAYERS

    "[OPTION] coordinates used to discretize the vertical"
    vertical_coordinates::VerticalCoordinates = SigmaCoordinates(; nlayers)

    # ARRAY TYPES (horizontal dimension in grid/spectral is flattened to 1D)
    SpectralVariable2D::Type{<:AbstractArray} = LowerTriangularArray{Complex{NF}, 1, ArrayType{Complex{NF}, 1}}
    SpectralVariable3D::Type{<:AbstractArray} = LowerTriangularArray{Complex{NF}, 2, ArrayType{Complex{NF}, 2}}
    SpectralVariable4D::Type{<:AbstractArray} = LowerTriangularArray{Complex{NF}, 3, ArrayType{Complex{NF}, 3}}
    GridVariable2D::Type{<:AbstractArray} = RingGrids.nonparametric_type(Grid){NF, 1, ArrayType{NF, 1}}
    GridVariable3D::Type{<:AbstractArray} = RingGrids.nonparametric_type(Grid){NF, 2, ArrayType{NF, 2}}
    GridVariable4D::Type{<:AbstractArray} = RingGrids.nonparametric_type(Grid){NF, 3, ArrayType{NF, 3}}
    ParticleVector::Type{<:AbstractArray} = ArrayType{Particle{NF}, 1}
end

function Base.show(io::IO, SG::SpectralGrid)
    (; NF, trunc, Grid, radius, nlat, npoints, nlayers, vertical_coordinates) = SG
    (; device, ArrayType) = SG
    (; nparticles) = SG

    # resolution information
    ave_resolution = sqrt(4π*radius^2/npoints)/1000  # in [km]
    s(x) = x > 1000 ? @sprintf("%i", x) : @sprintf("%.3g", x)

    println(io, "$(typeof(SG)):")
    println(io, "├ Spectral:   T$trunc LowerTriangularMatrix{Complex{$NF}}, radius = $radius m")
    println(io, "├ Grid:       $nlat-ring $Grid{$NF}, $npoints grid points")
    println(io, "├ Resolution: $(s(average_resolution))km (average)")
    if nparticles > 0
    println(io, "├ Particles:  $nparticles")
    end
    println(io, "├ Vertical:   $nlayers-layer $(typeof(vertical_coordinates))")
      print(io, "└ Device:     $(typeof(device)) using $ArrayType")
end

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct pulling in parameters from a SpectralGrid struct."""
function SpeedyTransforms.SpectralTransform(spectral_grid::SpectralGrid;
                                            one_more_degree::Bool = true,
                                            kwargs...)
    (; NF, Grid, trunc, nlat_half, nlayers, ArrayType) = spectral_grid
    return SpectralTransform(NF, trunc+one_more_degree, trunc, nlat_half; Grid, ArrayType, nlayers, kwargs...)
end