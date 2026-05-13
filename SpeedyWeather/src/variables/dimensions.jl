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
struct Grid3D <: AbstractVariableDim end
allocate(::AbstractVariable{Grid3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model))

"""Dimension for 3D land surface variables."""
struct Land3D <: AbstractVariableDim end
allocate(::AbstractVariable{Land3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.land))

"""Dimension for 3D ocean variables."""
struct Ocean3D <: AbstractVariableDim end
allocate(::AbstractVariable{Ocean3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model.ocean))

"""Dimension for 2D spectral variables."""
struct Spectral2D <: AbstractVariableDim end
allocate(::AbstractVariable{Spectral2D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable2D, model.spectral_grid.spectrum)

"""Dimension for 3D spectral variables with e.g. `n` vertical levels."""
@kwdef struct Spectral3D <: AbstractVariableDim
    n::Int = 0
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
allocate(v::AbstractVariable{Grid4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, get_nlayers(model), v.dims.n)

"""Dimension for 4D spectral variables with customizable extra dimension."""
@kwdef struct Spectral4D <: AbstractVariableDim
    n::Int = 1
end
allocate(v::AbstractVariable{Spectral4D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, get_nlayers(model), v.dims.n)

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

# Variable fusion support
# We may want to fuse a group of variables into a single parent variable to
# optimize performance by batching transforms.
#
# Two fuse families are supported:
#   :grid     — mix of Grid2D and Grid3D; parent is a Grid3D buffer
#   :spectral — mix of Spectral2D and Spectral3D; parent is a Spectral3D buffer
#
# Within a family, members can mix 2D and 3D freely:
#   - 3D members contribute `nlayers` slots and get a range-indexed view (keeps the layer dim)
#   - 2D members contribute 1 slot and get a scalar-indexed view (collapses the layer dim,
#     preserving the original 2D field/LTA semantics)
#
# `fused_slots` says how many parent-axis columns a member needs.
# `fuse_family` says which parent type the member is compatible with; fusion across families
# (e.g. mixing a Grid3D with a Spectral2D) is rejected in `build_fuse_parents`.

fused_slots(::Grid3D, model::AbstractModel) = get_nlayers(model)
fused_slots(::Grid2D, ::AbstractModel) = 1
fused_slots(::Spectral2D, ::AbstractModel) = 1
fused_slots(d::Spectral3D, model::AbstractModel) = d.n == 0 ? get_nlayers(model) : d.n

fuse_family(::Grid2D) = :grid
fuse_family(::Grid3D) = :grid
fuse_family(::Spectral2D) = :spectral
fuse_family(::Spectral3D) = :spectral
fuse_family(d::AbstractVariableDim) = error(
    "Fusion is not supported for dim type $(typeof(d)). " *
    "Supported dim types: Grid2D, Grid3D, Spectral2D, Spectral3D."
)

# Allocate the fused parent buffer for a group of variables sharing a fuse family.
# Variables are passed in declaration order; returns (parent, views, slots) where
# `views[i]` is the field/lta_view for `vars[i]` and `slots[i]` is its UnitRange{Int}
# along the fused axis of the parent. 2D members get scalar-indexed views (collapsed
# layer dim); 3D members get range-indexed views.
function allocate_fused(vars::AbstractVector{<:AbstractVariable}, model::AbstractModel)
    family = fuse_family(first(vars).dims)
    total = sum(fused_slots(v.dims, model) for v in vars)
    if family === :grid
        parent = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, total)
    else # :spectral (already validated upstream)
        parent = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, total)
    end
    views, slots = _split_views(parent, vars, model)

    return parent, views, slots
end

# Build views & slot ranges for a grid-family fuse group. Grid2D members use a scalar
# layer index so the resulting Field stays 2D (no layer axis); Grid3D members use a range.
function _split_views(parent, vars, model)
    views = Vector{Any}(undef, length(vars))
    slots = Vector{UnitRange{Int}}(undef, length(vars))
    offset = 0
    for (i, v) in enumerate(vars)
        n = fused_slots(v.dims, model)
        slots[i] = (offset + 1):(offset + n)
        views[i] = typeof(v.dims) <: Union{Grid2D, Spectral2D} ?
            wrapped_view(parent, :, offset + 1) :    # scalar → 2D array (no layer dim)
            wrapped_view(parent, :, slots[i])        # range  → 3D array (with layer dim)
        offset += n
    end
    @assert offset == size(parent.data, 2) "Fused grid parent has $offset slots assigned but parent has $(size(parent.data, 2)) layers"
    return views, slots
end
