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
const DEFAULT_NLAYERS_SOIL = 2

export SpectralGrid

"""
Defines the horizontal spectral resolution and corresponding grid and the
vertical coordinate for SpeedyWeather.jl. Options are
$(TYPEDFIELDS)

`nlat_half` and `npoints` should not be chosen but are derived from `trunc`,
`Grid` and `dealiasing`."""
struct SpectralGrid{
    SP,          # <: AbstractSpectrum
} <: AbstractSpectralGrid

    "[OPTION] number format used throughout the model"
    NF::Type{<:AbstractFloat}

    "[OPTION] device archictecture to run on"
    device::AbstractDevice

    "[OPTION] array type to use for all variables"
    ArrayType::Type{<:AbstractArray}

    "[DERIVED] Type of vector"
    VectorType::Type{<:AbstractVector}

    "[DERIVED] Type of matrix"
    MatrixType::Type{<:AbstractMatrix}

    "[DERIVED] Type of 3D array"
    TensorType::Type{<:AbstractArray}
    
    # HORIZONTAL SPECTRAL
    "[OPTION] horizontal resolution as the maximum degree of spherical harmonics"
    trunc::Int

    "[DERIVED] spectral space"
    spectrum::SP

    "[DERIVED] Type of spectral variable in 2D (horizontal only, flattened into 1D vector)"
    SpectralVariable2D::Type{<:AbstractArray}

    "[DERIVED] Type of spectral variable in 3D (horizontal only + e.g vertical, flattened into 2D matrix)"
    SpectralVariable3D::Type{<:AbstractArray}

    "[DERIVED] Type of spectral variable in 4D (horizontal only + e.g. vertical and time, flattened into 3D array)"
    SpectralVariable4D::Type{<:AbstractArray}
    
    # HORIZONTAL GRID
    "[OPTION] horizontal grid used for calculations in grid-point space"
    Grid::Type{<:AbstractGrid}

    "[DERIVED] Type of grid variable in 2D (horizontal only, flattened into 1D vector)"
    GridVariable2D::Type{<:AbstractArray}
    
    "[DERIVED] Type of grid variable in 3D (horizontal + e.g. vertical, flattened into 2D matrix)"
    GridVariable3D::Type{<:AbstractArray}
    
    "[DERIVED] Type of grid variable in 4D (horizontal + e.g. vertical + time, flattened into 3D array)"
    GridVariable4D::Type{<:AbstractArray}

    "[OPTION] how to match spectral with grid resolution: dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64

    # TODO move to planet?
    "[OPTION] radius of the sphere [m]"
    radius::Float64

    # PARTICLES
    "[OPTION] number of particles for particle advection [1]"
    nparticles::Int

    "[DERIVED] ArrayType of particle vector"
    ParticleVector::Type{<:AbstractArray}

    # SIZE OF GRID from trunc, Grid, dealiasing:
    "[DERIVED] number of latitude rings on one hemisphere (Equator incl)"
    nlat_half::Int

    "[DERIVED] number of latitude rings on both hemispheres"
    nlat::Int

    "[DERIVED] total number of grid points in the horizontal"
    npoints::Int

    # VERTICAL
    "[OPTION] number of vertical layers in the atmosphere"
    nlayers::Int

    "[OPTION] number of vertical layers in the soil/land"
    nlayers_soil::Int
end

function Base.show(io::IO, SG::SpectralGrid)
    (; NF, trunc, Grid, radius, nlat, npoints, nlayers, nlayers_soil) = SG
    (; device, ArrayType) = SG
    (; nparticles) = SG

    # resolution information
    average_resolution = sqrt(4π*radius^2/npoints)/1000  # in [km]
    s(x) = x > 1000 ? @sprintf("%i", x) : @sprintf("%.3g", x)

    println(io, "$(typeof(SG)):")
    println(io, "├ Spectral:   T$trunc LowerTriangularMatrix{Complex{$NF}}, radius = $radius m")
    println(io, "├ Grid:       $nlat-ring $Grid{$NF}, $npoints grid points")
    println(io, "├ Resolution: $(s(average_resolution))km (average)")
    nparticles > 0 &&
    println(io, "├ Particles:  $nparticles")
    println(io, "├ Vertical:   $nlayers-layer atmosphere, $nlayers_soil-layer land")
      print(io, "└ Device:     $(typeof(device)) using $ArrayType")
end

# Constructor that takes all [OPTION] parameters as keyword arguments
# and calculates all derived fields
function SpectralGrid(;
    NF::Type{<:AbstractFloat} = DEFAULT_NF,
    device::AbstractDevice = DEFAULT_DEVICE,
    ArrayType::Type{<:AbstractArray} = default_array_type(device),
    trunc::Int = DEFAULT_TRUNC,
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID,
    dealiasing::Real = 2.0,
    radius::Real = DEFAULT_RADIUS,
    nparticles::Int = 0,
    nlayers::Int = DEFAULT_NLAYERS,
    nlayers_soil::Int = DEFAULT_NLAYERS_SOIL
)
    # Convert numeric parameters to Float64
    dealiasing_f64 = Float64(dealiasing)
    radius_f64 = Float64(radius)
    
    # Calculate derived fields
    VectorType = ArrayType{NF, 1}
    MatrixType = ArrayType{NF, 2}
    TensorType = ArrayType{NF, 3}
    
    # Spectral space
    spectrum = Spectrum(trunc+2, trunc+1)
    
    # Spectral variable types
    SpectralVariable2D = LowerTriangularArray{Complex{NF}, 1, ArrayType{Complex{NF}, 1}, typeof(spectrum)}
    SpectralVariable3D = LowerTriangularArray{Complex{NF}, 2, ArrayType{Complex{NF}, 2}, typeof(spectrum)}
    SpectralVariable4D = LowerTriangularArray{Complex{NF}, 3, ArrayType{Complex{NF}, 3}, typeof(spectrum)}
    
    # Grid variable types
    GridVariable2D = RingGrids.nonparametric_type(Grid){NF, 1, ArrayType{NF, 1}}
    GridVariable3D = RingGrids.nonparametric_type(Grid){NF, 2, ArrayType{NF, 2}}
    GridVariable4D = RingGrids.nonparametric_type(Grid){NF, 3, ArrayType{NF, 3}}
    
    # Particle vector type
    ParticleVector = ArrayType{Particle{NF}, 1}
    
    # Grid size calculations
    nlat_half = SpeedyTransforms.get_nlat_half(trunc, dealiasing_f64)
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    npoints = RingGrids.get_npoints(Grid, nlat_half)
    
    # Create the SpectralGrid with all fields
    return SpectralGrid{typeof(spectrum)}(
        NF,
        device,
        ArrayType,
        VectorType,
        MatrixType,
        TensorType,
        trunc,
        spectrum,
        SpectralVariable2D,
        SpectralVariable3D,
        SpectralVariable4D,
        Grid,
        GridVariable2D,
        GridVariable3D,
        GridVariable4D,
        dealiasing_f64,
        radius_f64,
        nparticles,
        ParticleVector,
        nlat_half,
        nlat,
        npoints,
        nlayers,
        nlayers_soil
    )
end

# also allow spectral grid to be passed on as first an only positional argument to model constructors
(M::Type{<:AbstractModel})(SG::SpectralGrid; kwargs...) = M(spectral_grid=SG; kwargs...)

"""$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct pulling in parameters from a SpectralGrid struct."""
function SpeedyTransforms.SpectralTransform(spectral_grid::SpectralGrid;
                                            one_more_degree::Bool=true,
                                            kwargs...)
    (; NF, Grid, spectrum, nlat_half, nlayers, ArrayType) = spectral_grid
    if one_more_degree == false 
        return SpectralTransform(NF, Spectrum(spectrum.lmax-1, spectrum.mmax), nlat_half; Grid, ArrayType, nlayers, kwargs...)
    else 
        return SpectralTransform(NF, spectrum, nlat_half; Grid, ArrayType, nlayers, kwargs...)
    end 
end