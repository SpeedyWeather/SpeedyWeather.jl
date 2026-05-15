"""Dimension for the simulation clock and time tracking."""
struct ClockDim <: AbstractVariableDim end
allocate(::AbstractVariable{ClockDim}, model::AbstractModel) = Clock(model.architecture)

"""Dimension for scalar variables holding single numerical values as RefValue.
Default value is 0, but can be set via ScalarDim, e.g. ScalarDim(1) for a default value of 1."""
@kwdef struct ScalarDim{T} <: AbstractVariableDim
    value::T = zero(DEFAULT_NF)                                 # default value via ScalarDim(1)
end
allocate(v::AbstractVariable{<:ScalarDim}, model::AbstractModel) = Ref(convert(model.spectral_grid.NF, v.dims.value))

"""Dimension for 2D grid variables on the horizontal grid."""
struct Grid2D <: AbstractVariableDim end
allocate(::AbstractVariable{Grid2D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable2D, model.spectral_grid.grid)

"""Dimension for 3D grid variables on the horizontal grid with vertical levels."""
@kwdef struct Grid3D <: AbstractVariableDim
    n::Int = 0                                                  # length of 3rd dimension, use nlayers for 0
end
allocate(v::AbstractVariable{Grid3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n == 0 ? get_nlayers(model) : v.dims.n)

"""Dimension for 3D land surface variables."""
@kwdef struct Land3D <: AbstractVariableDim
    n::Int = 0                                                  # length of 3rd dimension, use nlayers for 0
end
allocate(v::AbstractVariable{Land3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n == 0 ? get_nlayers(model.land) : v.dims.n)

"""Dimension for 3D ocean variables."""
@kwdef struct Ocean3D <: AbstractVariableDim
    n::Int = 0                                                  # length of 3rd dimension, use nlayers for 0
end
allocate(v::AbstractVariable{Ocean3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n == 0 ? get_nlayers(model.ocean) : v.dims.n)

"""Dimension for 2D spectral variables."""
struct Spectral2D <: AbstractVariableDim end
allocate(::AbstractVariable{Spectral2D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable2D, model.spectral_grid.spectrum)

"""Dimension for 3D spectral variables with e.g. `n` vertical levels."""
@kwdef struct Spectral3D <: AbstractVariableDim
    n::Int = 0                                                  # length of 3rd dimension, use nlayers for 0
end
allocate(v::AbstractVariable{Spectral3D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, v.dims.n == 0 ? get_nlayers(model) : v.dims.n)

"""Dimension for 1D variables as a function of latitude only."""
struct Latitude1D <: AbstractVariableDim end
allocate(::AbstractVariable{Latitude1D}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, model.spectral_grid.nlat), 0)

"""Dimension for 1D vertical level variables."""
struct Vertical1D <: AbstractVariableDim end
allocate(::AbstractVariable{Vertical1D}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, get_nlayers(model)), 0)

"""Dimension for 4D grid variables with customizable extra dimension."""
@kwdef struct Grid4D <: AbstractVariableDim
    n::Int = 1                                                  # length of 4th dimension
end
allocate(v::AbstractVariable{Grid4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable4D, model.spectral_grid.grid, get_nlayers(model), v.dims.n)

"""Dimension for 4D spectral variables with customizable extra dimension."""
@kwdef struct Spectral4D <: AbstractVariableDim
    n::Int = 1                                                  # length of 4th dimension
end
allocate(v::AbstractVariable{Spectral4D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable4D, model.spectral_grid.spectrum, get_nlayers(model), v.dims.n)

"""Dimension for generic vector data of arbitrary length."""
@kwdef struct VectorDim <: AbstractVariableDim
    n::Int = 1
end
allocate(v::AbstractVariable{VectorDim}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, v.dims.n), 0)

"""Dimension for particle system vectors."""
@kwdef struct ParticleVectorDim <: AbstractVariableDim
    n::Int = 1
end
allocate(v::AbstractVariable{ParticleVectorDim}, model::AbstractModel) = zeros(model.spectral_grid.ParticleVectorType, v.dims.n)

"""Dimension for matrix data of arbitrary dimensions."""
@kwdef struct MatrixDim <: AbstractVariableDim
    m::Int = 1
    n::Int = 1
end
allocate(v::AbstractVariable{MatrixDim}, model::AbstractModel) = fill!(model.spectral_grid.MatrixType(undef, v.dims.m, v.dims.n), 0)

"""Dimension for spectral transform scratch memory."""
struct TransformScratchMemory <: AbstractVariableDim end
allocate(::AbstractVariable{TransformScratchMemory}, model::AbstractModel) = model.spectral_transform.scratch_memory

"""Dimension for particle locator tracking in space."""
@kwdef struct LocatorDim <: AbstractVariableDim
    n::Int = 1                                                  # number of locations to track, e.g. for particle advection
end
allocate(::AbstractVariable{LocatorDim}, model::AbstractModel) = RingGrids.AnvilLocator(model.spectral_grid.NF, model.particle_advection.nparticles; architecture = model.spectral_grid.architecture)
