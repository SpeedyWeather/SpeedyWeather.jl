abstract type AbstractSpectralGrid end

# computing
const DEFAULT_NF = Float32
const DEFAULT_ARCHITECTURE = CPU()
const DEFAULT_ARRAYTYPE = Array

# numerics
const DEFAULT_GRID = OctahedralGaussianGrid
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
    ArchitectureType,      # <: AbstractArchitecture
    SpectrumType,          # <: AbstractSpectrum
    GridType,              # <: AbstractGrid
} <: AbstractSpectralGrid

    "[OPTION] number format used throughout the model"
    NF::Type{<:AbstractFloat}

    "[OPTION] device architecture to run on"
    architecture::ArchitectureType 

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
    spectrum::SpectrumType

    "[DERIVED] Type of spectral variable in 2D (horizontal only, flattened into 1D vector)"
    SpectralVariable2D::Type{<:AbstractArray}

    "[DERIVED] Type of spectral variable in 3D (horizontal only + e.g vertical, flattened into 2D matrix)"
    SpectralVariable3D::Type{<:AbstractArray}

    "[DERIVED] Type of spectral variable in 4D (horizontal only + e.g. vertical and time, flattened into 3D array)"
    SpectralVariable4D::Type{<:AbstractArray}
    
    # SIZE OF GRID from trunc, Grid, dealiasing:
    "[OPTION] how to match spectral with grid resolution: dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64
    
    "[DERIVED] number of latitude rings on one hemisphere (Equator incl)"
    nlat_half::Int

    "[DERIVED] number of latitude rings on both hemispheres"
    nlat::Int

    "[DERIVED] total number of grid points in the horizontal"
    npoints::Int

    "[DERIVED] instance of horizontal grid used for calculations in grid-point space"
    grid::GridType

    "[OPTION] type of horizontal grid used for calculations in grid-point space"
    Grid::Type{<:AbstractGrid}

    "[DERIVED] Type of grid variable in 2D (horizontal only, flattened into 1D vector)"
    GridVariable2D::Type{<:AbstractArray}
    
    "[DERIVED] Type of grid variable in 3D (horizontal + e.g. vertical, flattened into 2D matrix)"
    GridVariable3D::Type{<:AbstractArray}

    "[DERIVED] Type of grid variable in 4D (horizontal + e.g. vertical + time, flattened into 3D array)"
    GridVariable4D::Type{<:AbstractArray}

    # PARTICLES
    "[OPTION] number of particles for particle advection [1]"
    nparticles::Int

    "[DERIVED] ArrayType of particle vector"
    ParticleVector::Type{<:AbstractArray}

    # VERTICAL
    "[OPTION] number of vertical layers in the atmosphere"
    nlayers::Int

    "[OPTION] number of vertical layers in the soil/land"
    nlayers_soil::Int
end

function Base.show(io::IO, SG::SpectralGrid)
    (; NF, trunc, grid, nlat, npoints, nlayers, nlayers_soil) = SG
    (; architecture, ArrayType) = SG
    (; nparticles) = SG
    Grid = nonparametric_type(grid)

    # resolution information
    radius = DEFAULT_RADIUS
    average_resolution = sqrt(4π*radius^2/npoints)/1000  # in [km]
    s(x) = x > 1000 ? @sprintf("%i", x) : @sprintf("%.3g", x)
    radius_str = @sprintf("%.0f", radius/1000)
    average_degrees = 360/sqrt(npoints*π)

    println(io, "SpectralGrid{Spectrum{...}, $Grid{...}}")
    println(io, "├ Number format: $NF")
    println(io, "├ Spectral:      T$trunc LowerTriangularMatrix")
    println(io, "├ Grid:          $nlat-ring $Grid, $npoints grid points")
    println(io, "├ Resolution:    $(s(average_degrees))°, $(s(average_resolution))km (at $(radius_str)km radius)")
    nparticles > 0 &&
    println(io, "├ Particles:     $nparticles")
    println(io, "├ Vertical:      $nlayers-layer atmosphere, $nlayers_soil-layer land")
    print(io,   "└ Architecture:  $architecture using $ArrayType")
end

# Constructor that takes all [OPTION] parameters as keyword arguments
# and calculates all derived fields
"""$(TYPEDSIGNATURES) 
Initialize a SpectralGrid from a given truncation and all [OPTION] parameters of SpectralGrid."""
function SpectralGrid(;
    NF::Type{<:AbstractFloat} = DEFAULT_NF,
    architecture::Union{AbstractArchitecture, Type{<:AbstractArchitecture}} = DEFAULT_ARCHITECTURE,
    trunc::Int = DEFAULT_TRUNC,
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID,
    dealiasing::Real = 2,
    nparticles::Int = 0,
    nlayers::Int = DEFAULT_NLAYERS,
    nlayers_soil::Int = DEFAULT_NLAYERS_SOIL
)

    # Convert architecture to instance if it is a type
    if architecture isa Type
        architecture = architecture()
    end

    # grid
    nlat_half = SpeedyTransforms.get_nlat_half(trunc, dealiasing)
    grid = Grid(nlat_half, architecture)

    # default dealiasing or user-defined one? 
    dealiasing = SpeedyTransforms.get_dealiasing(trunc, grid.nlat_half)
    
    # Spectral space
    spectrum = Spectrum(trunc+2, trunc+1, architecture=architecture)

    # Create the SpectralGrid with all fields
    return SpectralGrid(NF, spectrum, grid, dealiasing, nparticles, nlayers, nlayers_soil)
end

"""
$(TYPEDSIGNATURES)
Initialize a SpectralGrid from a given grid.
"""
function SpectralGrid(grid::AbstractGrid,
    NF::Type{<:AbstractFloat} = DEFAULT_NF,
    dealiasing::Real = 2,
    nparticles::Int = 0,
    nlayers::Int = DEFAULT_NLAYERS,
    nlayers_soil::Int = DEFAULT_NLAYERS_SOIL
)
    architecture = grid.architecture
    
    trunc = get_truncation(grid, dealiasing)
    nlat_half = get_nlat_half(grid)
    nlat = RingGrids.get_nlat(grid)
    npoints = RingGrids.get_npoints(grid)

    spectrum = Spectrum(trunc+2, trunc+1, architecture=architecture)

    return SpectralGrid(NF, spectrum, grid, dealiasing, nparticles, nlayers, nlayers_soil)
end

# low level constructor, not intended to be used directly by users
function SpectralGrid(NF::Type{<:AbstractFloat},
    spectrum::Spectrum, 
    grid::AbstractGrid, 
    dealiasing,
    nparticles,
    nlayers::Int,
    nlayers_soil::Int
    )
    @assert spectrum.architecture == grid.architecture "Architecture of grid and spectrum must match"

    architecture = spectrum.architecture

    # grid
    nlat_half = SpeedyTransforms.get_nlat_half(grid)
    nlat = RingGrids.get_nlat(grid)
    npoints = RingGrids.get_npoints(grid)

    # Convert numeric parameters to Float64
    dealiasing_f64 = Float64(dealiasing)

    # Calculate derived fields
    ArrayType = array_type(architecture)
    VectorType = array_type(architecture, NF, 1)
    MatrixType = array_type(architecture, NF, 2)
    TensorType = array_type(architecture, NF, 3)
    
    # Spectral variable types
    SpectralVariable2D = LowerTriangularArray{Complex{NF}, 1, array_type(architecture, Complex{NF}, 1), typeof(spectrum)}
    SpectralVariable3D = LowerTriangularArray{Complex{NF}, 2, array_type(architecture, Complex{NF}, 2), typeof(spectrum)}
    SpectralVariable4D = LowerTriangularArray{Complex{NF}, 3, array_type(architecture, Complex{NF}, 3), typeof(spectrum)}

    # Grid variable types
    GridVariable2D = Field{NF, 1, array_type(architecture, NF, 1), typeof(grid)}
    GridVariable3D = Field{NF, 2, array_type(architecture, NF, 2), typeof(grid)}
    GridVariable4D = Field{NF, 3, array_type(architecture, NF, 3), typeof(grid)}
    
    # Particle vector type
    ParticleVector = array_type(architecture, Particle{NF}, 1)
    
    # Create the SpectralGrid with all fields
    return SpectralGrid{typeof(architecture), typeof(spectrum), typeof(grid)}(
        NF,
        architecture,
        ArrayType,
        VectorType,
        MatrixType,
        TensorType,
        truncation(spectrum),
        spectrum,
        SpectralVariable2D,
        SpectralVariable3D,
        SpectralVariable4D,
        dealiasing_f64,
        nlat_half,
        nlat,
        npoints,
        grid,
        nonparametric_type(grid),
        GridVariable2D,
        GridVariable3D,
        GridVariable4D,
        nparticles,
        ParticleVector,
        nlayers,
        nlayers_soil
    )
end

# also allow spectral grid to be passed on as first and only positional argument to model constructors
(M::Type{<:AbstractModel})(SG::SpectralGrid; kwargs...) = M(; spectral_grid=SG, kwargs...)

"""$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct pulling in parameters from a SpectralGrid struct."""
function SpeedyTransforms.SpectralTransform(spectral_grid::SpectralGrid;
                                            one_more_degree::Bool=true,
                                            kwargs...)
    (; NF, spectrum, grid, nlayers, ArrayType) = spectral_grid
    (; lmax, mmax, architecture) = spectrum
    spectrum = one_more_degree == false ? Spectrum(lmax-1, mmax; architecture) : spectrum
    return SpectralTransform(spectrum, grid; NF, ArrayType, nlayers, kwargs...)
end

# because model components can be `nothing`, their constructor being `Nothing()`
# we also allow `::SpectralGrid` as the first argument
Base.Nothing(::SpectralGrid) = Nothing()