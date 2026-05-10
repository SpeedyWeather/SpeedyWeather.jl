"""Dimension for the simulation clock and time tracking."""
struct ClockDim <: AbstractVariableDim end
Base.zero(::AbstractVariable{ClockDim}, model::AbstractModel) = Clock(model.architecture)

"""Dimension for scalar variables holding single numerical values as RefValue.
Default value is 0, but can be set via ScalarDim, e.g. ScalarDim(1) for a default value of 1."""
@kwdef struct ScalarDim{T} <: AbstractVariableDim
    value::T = zero(DEFAULT_NF)                                 # default value via ScalarDim(1)
end
Base.zero(v::AbstractVariable{<:ScalarDim}, model::AbstractModel) = Ref(convert(model.spectral_grid.NF, v.dims.value))

"""Dimension for 2D grid variables on the horizontal grid."""
struct Grid2D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Grid2D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable2D, model.spectral_grid.grid)

"""Dimension for 3D grid variables on the horizontal grid with vertical levels."""
struct Grid3D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Grid3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model))

"""Dimension for 3D land surface variables."""
struct Land3D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Land3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.land))

"""Dimension for 3D ocean variables."""
struct Ocean3D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Ocean3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.ocean))

"""Dimension for 2D spectral variables."""
struct Spectral2D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Spectral2D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable2D, model.spectral_grid.spectrum)

"""Dimension for 3D spectral variables with e.g. `n` vertical levels."""
@kwdef struct Spectral3D <: AbstractVariableDim
    n::Int = 0
end
Base.zero(v::AbstractVariable{Spectral3D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, v.dims.n == 0 ? get_nlayers(model) : v.dims.n)

"""Dimension for 1D variables as a function of latitude only."""
struct Latitude1D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Latitude1D}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, model.spectral_grid.nlat), 0)

"""Dimension for 1D vertical level variables."""
struct Vertical1D <: AbstractVariableDim end
Base.zero(::AbstractVariable{Vertical1D}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, get_nlayers(model)), 0)

"""Dimension for 4D grid variables with customizable extra dimension."""
@kwdef struct Grid4D <: AbstractVariableDim
    n::Int = 1                                                  # length of 4th dimension
end
Base.zero(v::AbstractVariable{Grid4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model), v.dims.n)

"""Dimension for 4D spectral variables with customizable extra dimension."""
@kwdef struct Spectral4D <: AbstractVariableDim
    n::Int = 1
end
Base.zero(v::AbstractVariable{Spectral4D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, get_nlayers(model), v.dims.n)

"""Dimension for generic vector data of arbitrary length."""
@kwdef struct VectorDim <: AbstractVariableDim
    n::Int = 1
end
Base.zero(v::AbstractVariable{VectorDim}, model::AbstractModel) = fill!(model.spectral_grid.VectorType(undef, v.dims.n), 0)

"""Dimension for particle system vectors."""
@kwdef struct ParticleVectorDim <: AbstractVariableDim
    n::Int = 1
end
Base.zero(v::AbstractVariable{ParticleVectorDim}, model::AbstractModel) = zeros(model.spectral_grid.ParticleVectorType, v.dims.n)

"""Dimension for matrix data of arbitrary dimensions."""
@kwdef struct MatrixDim <: AbstractVariableDim
    m::Int = 1
    n::Int = 1
end
Base.zero(v::AbstractVariable{MatrixDim}, model::AbstractModel) = fill!(model.spectral_grid.MatrixType(undef, v.dims.m, v.dims.n), 0)

"""Dimension for spectral transform scratch memory."""
struct TransformScratchMemory <: AbstractVariableDim end
Base.zero(::AbstractVariable{TransformScratchMemory}, model::AbstractModel) = model.spectral_transform.scratch_memory

"""Dimension for particle locator tracking in space."""
@kwdef struct LocatorDim <: AbstractVariableDim
    n::Int = 1                                                  # number of locations to track, e.g. for particle advection
end
Base.zero(::AbstractVariable{LocatorDim}, model::AbstractModel) = RingGrids.AnvilLocator(model.spectral_grid.NF, model.particle_advection.nparticles; architecture = model.spectral_grid.architecture)

# Variable fusion support
# We may want to fuse a group of variables into a single parent variable to 
# optimize performance by batching transforms. 
# Here we define how many "slots" along the layer axis a variable contributes
# to a fused parent, and how to allocate the parent + per-member views.
# Currently supports Grid3D, Grid2D, Spectral3D, Spectral2D. Variables in a
# fused group must share the same dim *type* (validated in allocate).

fused_slots(::Grid3D, model::AbstractModel) = get_nlayers(model)
fused_slots(::Grid2D, ::AbstractModel) = 1
fused_slots(::Spectral2D, ::AbstractModel) = 1
fused_slots(d::Spectral3D, model::AbstractModel) = d.n == 0 ? get_nlayers(model) : d.n

# allocate the fused parent buffer for a group of variables sharing a dim type.
# variables are passed in declaration order; the function returns (parent, views)
# where `views` is a vector of per-member field/lta_view aligned with `vars`.
function allocate_fused(vars::AbstractVector{<:AbstractVariable{Grid3D}}, model::AbstractModel)
    total = sum(fused_slots(v.dims, model) for v in vars)
    parent = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, total)
    return parent, _split_views_grid3d(parent, vars, model)
end

function allocate_fused(vars::AbstractVector{<:AbstractVariable{Grid2D}}, model::AbstractModel)
    # each Grid2D contributes 1 layer; parent is a Grid3D with N layers
    parent = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, length(vars))
    views = [field_view(parent, :, k) for k in 1:length(vars)]
    return parent, views
end

function allocate_fused(vars::AbstractVector{<:AbstractVariable{Spectral3D}}, model::AbstractModel)
    total = sum(fused_slots(v.dims, model) for v in vars)
    parent = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, total)
    return parent, _split_views_spectral3d(parent, vars, model)
end

function allocate_fused(vars::AbstractVector{<:AbstractVariable{Spectral2D}}, model::AbstractModel)
    parent = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, length(vars))
    views = [lta_view(parent, :, k) for k in 1:length(vars)]
    return parent, views
end

# fallback: clear error message for unsupported dim types
allocate_fused(::AbstractVector{<:AbstractVariable{D}}, ::AbstractModel) where {D} =
    error("Fusion is not implemented for dim type $D. Supported: Grid2D, Grid3D, Spectral2D, Spectral3D.")

function _split_views_grid3d(parent, vars, model)
    views = Vector{Any}(undef, length(vars))
    offset = 0
    for (i, v) in enumerate(vars)
        n = fused_slots(v.dims, model)
        views[i] = field_view(parent, :, (offset + 1):(offset + n))
        offset += n
    end
    return views
end

function _split_views_spectral3d(parent, vars, model)
    views = Vector{Any}(undef, length(vars))
    offset = 0
    for (i, v) in enumerate(vars)
        n = fused_slots(v.dims, model)
        views[i] = lta_view(parent, :, (offset + 1):(offset + n))
        offset += n
    end
    return views
end

