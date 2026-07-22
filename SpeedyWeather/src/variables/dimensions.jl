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
struct Grid2D <: AbstractVariableDim2D end
allocate(::AbstractVariable{Grid2D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable2D, model.spectral_grid.grid)
const GridXY = Grid2D

"""Dimension for 3D grid variables on the horizontal grid with an unspecified 3rd dimension
of length `n` (default 1, a trailing singleton dimension). Carries a neutral tag; prefer
`GridXYZ` or `GridXYT` when the 3rd dimension is specifically vertical or time."""
@kwdef struct Grid3D <: AbstractVariableDim3D
    n::Int = 1                                                  # length of 3rd dimension
end
allocate(v::AbstractVariable{Grid3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n)

"""Dimension for 3D grid variables on the horizontal grid with vertical levels (nlayers)."""
struct GridXYZ <: AbstractVariableDim3D end
allocate(::AbstractVariable{GridXYZ}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZ, model.spectral_grid.grid, get_nlayers(model))

"""Dimension for 3D grid variables on the horizontal grid with `n` time steps."""
@kwdef struct GridXYT <: AbstractVariableDim3D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{GridXYT}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYT, model.spectral_grid.grid, v.dims.n)

"""Dimension for 3D land surface variables, horizontal + vertical (nlayers)."""
struct LandXYZ <: AbstractVariableDim3D end
allocate(::AbstractVariable{LandXYZ}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZ, model.spectral_grid.grid, get_nlayers(model.land))
const Land3D = LandXYZ

"""Dimension for 3D land surface variables, horizontal + `n` time steps."""
@kwdef struct LandXYT <: AbstractVariableDim3D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{LandXYT}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYT, model.spectral_grid.grid, v.dims.n)

"""Dimension for 4D land surface variables, horizontal + vertical (nlayers) + `n` time steps."""
@kwdef struct Land4D <: AbstractVariableDim4D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{Land4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZT, model.spectral_grid.grid, get_nlayers(model.land), v.dims.n)
const LandXYZT = Land4D

"""Dimension for 3D ocean variables with an unspecified 3rd dimension of length `n`
(default 1, a trailing singleton dimension). Carries a neutral tag; prefer `OceanXYZ` or
`OceanXYT` when the 3rd dimension is specifically vertical or time."""
@kwdef struct Ocean3D <: AbstractVariableDim3D
    n::Int = 1                                                  # length of 3rd dimension
end
allocate(v::AbstractVariable{Ocean3D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable3D, model.spectral_grid.grid, v.dims.n)

"""Dimension for 3D ocean variables, horizontal + vertical (ocean layers)."""
struct OceanXYZ <: AbstractVariableDim3D end
allocate(::AbstractVariable{OceanXYZ}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZ, model.spectral_grid.grid, get_nlayers(model.ocean))

"""Dimension for 3D ocean variables, horizontal + `n` time steps."""
@kwdef struct OceanXYT <: AbstractVariableDim3D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{OceanXYT}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYT, model.spectral_grid.grid, v.dims.n)

"""Dimension for 4D ocean variables, horizontal + vertical (ocean layers) + `n` time steps."""
@kwdef struct Ocean4D <: AbstractVariableDim4D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{Ocean4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZT, model.spectral_grid.grid, get_nlayers(model.ocean), v.dims.n)
const OceanXYZT = Ocean4D

"""Dimension for 2D spectral variables."""
struct Spectral2D <: AbstractVariableDim2D end
allocate(::AbstractVariable{Spectral2D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable2D, model.spectral_grid.spectrum)

"""Dimension for 3D spectral variables with an unspecified 3rd dimension of length `n`
(default 1, a trailing singleton dimension). Carries a neutral tag; prefer `SpectralXYZ`
or `SpectralXYT` when the 3rd dimension is specifically vertical or time."""
@kwdef struct Spectral3D <: AbstractVariableDim3D
    n::Int = 1                                                  # length of 3rd dimension
end
allocate(v::AbstractVariable{Spectral3D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable3D, model.spectral_grid.spectrum, v.dims.n)

"""Dimension for 3D spectral variables, horizontal + vertical (nlayers)."""
struct SpectralXYZ <: AbstractVariableDim3D end
allocate(::AbstractVariable{SpectralXYZ}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariableXYZ, model.spectral_grid.spectrum, get_nlayers(model))

"""Dimension for 3D spectral variables, horizontal + `n` time steps."""
@kwdef struct SpectralXYT <: AbstractVariableDim3D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{SpectralXYT}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariableXYT, model.spectral_grid.spectrum, v.dims.n)

"""Dimension for 1D variables as a function of latitude only."""
struct Latitude1D <: AbstractVariableDim end
allocate(::AbstractVariable{Latitude1D}, model::AbstractModel) = on_architecture(model.spectral_grid.architecture, zeros(model.spectral_grid.NF, model.spectral_grid.nlat))

"""Dimension for 1D vertical level variables."""
struct Vertical1D <: AbstractVariableDim end
allocate(::AbstractVariable{Vertical1D}, model::AbstractModel) = on_architecture(model.spectral_grid.architecture, zeros(model.spectral_grid.NF, get_nlayers(model)))

"""Dimension for 4D grid variables with customizable extra dimension (neutral tag). Prefer
`GridXYZT` when the trailing dimensions are specifically vertical + time."""
@kwdef struct Grid4D <: AbstractVariableDim4D
    n::Int = 1                                                  # length of 4th dimension
end
allocate(v::AbstractVariable{Grid4D}, model::AbstractModel) = zeros(model.spectral_grid.GridVariable4D, model.spectral_grid.grid, get_nlayers(model), v.dims.n)

"""Dimension for 4D grid variables, horizontal + vertical (nlayers) + `n` time steps."""
@kwdef struct GridXYZT <: AbstractVariableDim4D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{GridXYZT}, model::AbstractModel) = zeros(model.spectral_grid.GridVariableXYZT, model.spectral_grid.grid, get_nlayers(model), v.dims.n)

"""Dimension for 4D spectral variables with customizable extra dimension (neutral tag). Prefer
`SpectralXYZT` when the trailing dimensions are specifically vertical + time."""
@kwdef struct Spectral4D <: AbstractVariableDim4D
    n::Int = 1                                                  # length of 4th dimension
end
allocate(v::AbstractVariable{Spectral4D}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariable4D, model.spectral_grid.spectrum, get_nlayers(model), v.dims.n)

"""Dimension for 4D spectral variables, horizontal + vertical (nlayers) + `n` time steps."""
@kwdef struct SpectralXYZT <: AbstractVariableDim4D
    n::Int = 1                                                  # length of time dimension
end
allocate(v::AbstractVariable{SpectralXYZT}, model::AbstractModel) = zeros(model.spectral_grid.SpectralVariableXYZT, model.spectral_grid.spectrum, get_nlayers(model), v.dims.n)

"""Dimension for generic vector data of arbitrary length."""
@kwdef struct VectorDim <: AbstractVariableDim
    n::Int = 1
end
allocate(v::AbstractVariable{VectorDim}, model::AbstractModel) = on_architecture(model.spectral_grid.architecture, zeros(model.spectral_grid.NF, v.dims.n))

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
allocate(v::AbstractVariable{MatrixDim}, model::AbstractModel) = on_architecture(model.spectral_grid.architecture, zeros(model.spectral_grid.NF, v.dims.m, v.dims.n))

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
#   :grid     — mix of Grid2D, Grid3D, Grid4D members
#   :spectral — mix of Spectral2D, Spectral3D, Spectral4D members
#
# Parent rank: 3D by default. As soon as any 4D member is in the group, the parent
# becomes 4D and inherits its trailing dim size `n` (timesteps/time step cache) from the 4D members.
#
# Concatenation happens along the *layer* axis (dim 2 of the parent's `.data`). The
# trainling/steps dim (when present) is shared by all members of the group, which preserves
# contiguity for per-step views: `view(parent.data, :, :, step)` is contiguous in memory.
#
# Member shapes (view rank) depending on parent rank:
#
#   parent 3D (no time steps)         | parent 4D `(npoints, total_slots, nsteps)`
#   ----------------------------------+---------------------------------------------
#   Grid2D     → 1 slot,  scalar idx  | DISALLOWED because no time step / n step
#                 view (npoints,)     |
#   Grid3D     → n slots, range       | 1 slot, scalar layer-idx, trailing `:`
#                 view (npoints, n)   → view (npoints, nsteps)
#   Grid4D     → forces parent to 4D  | nlayers slots, range, trailing `:`
#                                     → view (npoints, nlayers, nsteps)
#
# Same rules apply to the spectral family (Spectral2D/3D/4D) with `npoints` → `lm`.
#
# `fused_slots` says how many parent layer-axis columns a member needs. It is
# *parent-rank-dependent* via the `parent_is_4d` flag: a Grid3D member uses
# `n` slots in a 3D parent but only 1 slot in a 4D parent (because it is encoding
# a 2D field with steps for the solver in this case).
# `fuse_family` says which parent family the member is compatible with; fusion across
# families is rejected in `build_fuse_parents`.

fuse_family(::Grid2D) = :grid
fuse_family(::Grid3D) = :grid
fuse_family(::GridXYZ) = :grid
fuse_family(::GridXYT) = :grid
fuse_family(::Grid4D) = :grid
fuse_family(::GridXYZT) = :grid
fuse_family(::Spectral2D) = :spectral
fuse_family(::Spectral3D) = :spectral
fuse_family(::SpectralXYZ) = :spectral
fuse_family(::SpectralXYT) = :spectral
fuse_family(::Spectral4D) = :spectral
fuse_family(::SpectralXYZT) = :spectral
fuse_family(d::AbstractVariableDim) = error(
    "Fusion is not supported for dim type $(typeof(d)). " *
    "Supported dim types: Grid2D, Grid3D, GridXYZ, GridXYT, Grid4D, GridXYZT, " *
    "Spectral2D, Spectral3D, SpectralXYZ, SpectralXYT, Spectral4D, SpectralXYZT."
)

# Whether a member's dim forces the parent to be 4D.
is_fuse_4d(::Grid4D) = true
is_fuse_4d(::GridXYZT) = true
is_fuse_4d(::Spectral4D) = true
is_fuse_4d(::SpectralXYZT) = true
is_fuse_4d(::AbstractVariableDim) = false

# Whether a member's dim is a "horizontal-only" 2D field/LTA (no layer dim).
is_fuse_2d(::Grid2D) = true
is_fuse_2d(::Spectral2D) = true
is_fuse_2d(::AbstractVariableDim) = false

# Trailing-dim size for 4D members, how many n_steps
fuse_trailing_n(d::Grid4D) = d.n
fuse_trailing_n(d::GridXYZT) = d.n
fuse_trailing_n(d::Spectral4D) = d.n
fuse_trailing_n(d::SpectralXYZT) = d.n

# Layer-axis slot count for a member, depending on whether the parent is 3D or 4D.
# In a 4D parent every non-4D member collapses to 1 layer slot (it shares the trailing
# dim with the 4D members); only 4D members keep `nlayers` slots.
function fused_slots(d::AbstractVariableDim, model::AbstractModel; parent_is_4d::Bool = false)
    if parent_is_4d
        is_fuse_4d(d) ? get_nlayers(model) : 1
    else
        is_fuse_2d(d) ? 1 :
        d isa Spectral3D ? d.n :
        d isa SpectralXYT ? d.n :
        d isa Grid3D ? d.n :
        d isa GridXYT ? d.n :
        get_nlayers(model)
    end
end

# Allocate the fused parent buffer for a group of variables sharing a fuse family.
# Variables are passed in declaration order; returns (parent, views, slots) where
# `views[i]` is the field/lta_view for `vars[i]` and `slots[i]` is its UnitRange{Int}
# along the layer axis of the parent.
function allocate_fused(vars::AbstractVector{<:AbstractVariable}, model::AbstractModel)
    family = fuse_family(first(vars).dims)
    parent_is_4d, n_steps = _fuse_rank_and_n(vars)

    total = sum(fused_slots(v.dims, model; parent_is_4d) for v in vars)
    NF = model.spectral_grid.NF
    if family === :grid
        parent = parent_is_4d ?
            zeros(NF, model.spectral_grid.grid, ArrayDimensions.XYZT(), total, n_steps) :
            zeros(NF, model.spectral_grid.grid, ArrayDimensions.XYZ(), total)
    else
        parent = parent_is_4d ?
            zeros(Complex{NF}, model.spectral_grid.spectrum, ArrayDimensions.LMZT(), total, n_steps) :
            zeros(Complex{NF}, model.spectral_grid.spectrum, ArrayDimensions.LMZ(), total)
    end
    views, slots = _split_views(parent, vars, model, parent_is_4d)

    return parent, views, slots
end

# Decide the parent's rank for a fuse group and (if 4D) the shared trailing dim size.
# Enforces:
#   - all 4D members in the group share the same `n`
#   - if the group has any 4D member, no 2D member is allowed (2D members would need a
#     separate parent rank; they have to live in their own fuse group)
function _fuse_rank_and_n(vars::AbstractVector{<:AbstractVariable})
    fourD_members = filter(v -> is_fuse_4d(v.dims), vars)
    isempty(fourD_members) && return (false, 0)
    ns_set = unique(fuse_trailing_n(v.dims) for v in fourD_members)
    length(ns_set) == 1 || error(
        "Fuse group $(first(vars).fuse) has 4D members with mixed trailing-dim sizes $(collect(ns_set)). " *
        "All 4D members of a fuse group must share the same `n`."
    )
    for v in vars
        is_fuse_2d(v.dims) && error(
            "Fuse group: variable `$(v.name)` has dim $(typeof(v.dims)) but the group " *
            "contains 4D members which force a 4D parent. 2D members cannot be fused into " *
            "a 4D parent — give them their own fuse symbol."
        )
    end
    return (true, first(ns_set))
end

# Build views & slot ranges for a fuse group.
#   - 3D parent: 2D members get scalar layer-idx (drops layer dim); 3D members get a range.
#   - 4D parent: 3D members get scalar layer-idx + trailing `:` (so view is (npoints, n));
#                4D members get a range + trailing `:` (so view is (npoints, nlayers, n)).
function _split_views(parent, vars, model, parent_is_4d::Bool)
    views = Vector{Any}(undef, length(vars))
    slots = Vector{UnitRange{Int}}(undef, length(vars))
    offset = 0
    for (i, v) in enumerate(vars)
        n = fused_slots(v.dims, model; parent_is_4d)
        slots[i] = (offset + 1):(offset + n)
        views[i] = if parent_is_4d
            if is_fuse_4d(v.dims)
                wrapped_view(parent, :, slots[i], :)         # 4D member → (npoints, nlayers, n)
            else
                wrapped_view(parent, :, offset + 1, :)       # 3D-in-4D member → (npoints, n)
            end
        else
            is_fuse_2d(v.dims) ?
                wrapped_view(parent, :, offset + 1) :        # 2D → scalar, drop layer dim
                wrapped_view(parent, :, slots[i])            # 3D → range, keep layer dim
        end
        offset += n
    end
    @assert offset == size(parent.data, 2) "Fused parent has $offset slots assigned but axis 2 has size $(size(parent.data, 2))"
    return views, slots
end
