# as a fallback use one step for prognostic variables
prognostic_steps(::AbstractTimeStepper) = 1
prognostic_grid_steps(TS::AbstractTimeStepper) = prognostic_steps(TS)
prognostic_spectral_steps(TS::AbstractTimeStepper) = prognostic_steps(TS)

# as a fallback use one step for tendencies
tendency_steps(::AbstractTimeStepper) = 1
tendency_grid_steps(TS::AbstractTimeStepper) = tendency_steps(TS)
tendency_spectral_steps(TS::AbstractTimeStepper) = tendency_steps(TS)

# also allow for the model to be dispatched over
prognostic_steps(TS::AbstractTimeStepper, ::AbstractModel) = prognostic_steps(TS)
prognostic_grid_steps(TS::AbstractTimeStepper, ::AbstractModel) = prognostic_grid_steps(TS)
prognostic_spectral_steps(TS::AbstractTimeStepper, ::AbstractModel) = prognostic_spectral_steps(TS)

tendency_steps(TS::AbstractTimeStepper, ::AbstractModel) = tendency_steps(TS)
tendency_grid_steps(TS::AbstractTimeStepper, ::AbstractModel) = tendency_grid_steps(TS)
tendency_spectral_steps(TS::AbstractTimeStepper, ::AbstractModel) = tendency_spectral_steps(TS)

const DEFAULT_NSTEPS = (
    prognostic_grid = prognostic_grid_steps(DummyTimeStepper()),
    prognostic_spectral = prognostic_spectral_steps(DummyTimeStepper()),
    tendency_grid = tendency_grid_steps(DummyTimeStepper()),
    tendency_spectral = tendency_spectral_steps(DummyTimeStepper()),
)

function get_nsteps(time_stepping::AbstractTimeStepper, model::AbstractModel)
    return (;
        prognostic_grid = prognostic_grid_steps(time_stepping, model),
        prognostic_spectral = prognostic_spectral_steps(time_stepping, model),
        tendency_grid = tendency_grid_steps(time_stepping, model),
        tendency_spectral = tendency_spectral_steps(time_stepping, model),
    )
end

"""$(TYPEDSIGNATURES)
Get all steps of a variable as a tuple of views (wrapped into the same type as the input variable)
as defined by `get_step`. "Steps" refer to the step dimension, i.e. the time dimension `T` in the
variable's `ArrayDimensions` (e.g. `XYT`, `LMZT`), for prognostic variables e.g. used for the
leapfrog time step. Variables without a time dimension have no steps to splat, so a 1-tuple
holding the full variable (as a view, consistent with `get_step`) is returned."""
get_steps

# variables WITH a step dimension: splat the step dimension into a tuple, one view per step
# variables WITHOUT one: nothing to splat, a single "step" that is the full variable as a view
get_steps(var::Union{AbstractField, LowerTriangularArray}) = ntuple(step -> get_step(var, step), nsteps(var))

# plain arrays (e.g. unwrapped inside GPU kernels) carry no dimension information, so the
# last dimension is assumed to be the step dimension, as in `get_step`
get_steps(var::AbstractArray{T, 1}) where {T} = (var,)
get_steps(var::AbstractArray{T, 2}) where {T} = ntuple(step -> get_step(var, step), size(var, 2))
get_steps(var::AbstractArray{T, 3}) where {T} = ntuple(step -> get_step(var, step), size(var, 3))

"""$(TYPEDSIGNATURES)
Number of steps of a variable, i.e. the length of its step dimension (the time dimension `T`
in the variable's `ArrayDimensions`). Variables without a time dimension have a single step."""
# dispatch on the concrete dimension types that have a time dimension, NOT on the
# ...WithTime unions: their `Dims <: DimensionsWithTime` bound does not make the signature
# more specific than the unbounded one, so the two would be ambiguous and resolved by
# definition order rather than by the presence of a time dimension.
@inline nsteps(var::AbstractField{T, 1, A, G, ArrayDimensions.XYT}) where {T, A, G} = size(var, ndims(var))
@inline nsteps(var::AbstractField{T, 2, A, G, ArrayDimensions.XYT}) where {T, A, G} = size(var, ndims(var))
@inline nsteps(var::AbstractField{T, 3, A, G, ArrayDimensions.XYZT}) where {T, A, G} = size(var, ndims(var))
@inline nsteps(var::LowerTriangularArray{T, 2, A, S, ArrayDimensions.LMT}) where {T, A, S} = size(var, ndims(var))
@inline nsteps(var::LowerTriangularArray{T, 3, A, S, ArrayDimensions.LMZT}) where {T, A, S} = size(var, ndims(var))
@inline nsteps(::AbstractField) = 1
@inline nsteps(::LowerTriangularArray) = 1

"""$(TYPEDSIGNATURES)
Get the first `N` steps of a variable as a tuple of views with a COMPILE-TIME length,
fully unrolled via `Val(N)` (branchless, in contrast to the runtime-length methods above
whose `ntuple(f, ::Int)` returns a small union of tuple types; that union breaks Enzyme's
type analysis on Julia ≥ 1.11, see https://github.com/EnzymeAD/Enzyme.jl/issues/3275).
NOTE: in Enzyme-differentiated code that launches kernels on the views, even a compile-time
tuple of large step views (e.g. of a `LowerTriangularArray`) can exceed Enzyme's type
analysis — there, avoid the tuple altogether and bind steps individually via `get_step`."""
get_steps(var::AbstractArray, ::Val{N}) where {N} = ntuple(step -> get_step(var, step), Val(N))

export get_step

"""$(TYPEDSIGNATURES)
Select the step dimension from a variable, when no step is provided as 2nd argument select the
last step as this typically represents the "current" step (and not any previous ones). But this
depends on the time stepping a variable with step dimension was created for. Variables without a
step dimension have a single step, so `get_step(var)` returns the full variable as a view."""
get_step(var) = get_step(var, nsteps(var))

# plain arrays carry no dimension information so the last dimension is assumed to be the step one
get_step(var::AbstractArray) = get_step(var, size(var, ndims(var)))

# Plain Arrays
# Inside GPU kernels `Adapt.adapt_structure(to, field::AbstractField) = adapt(to, field.data)`
# unwraps a Field to its bare device array, so `get_step` must also work on plain arrays
# (otherwise the device-side MethodError aborts kernel compilation). Without the dimension
# information of a Field/LowerTriangularArray, `step` here always indexes the last dimension;
# for arrays without an explicit step dimension `step` must be 1 (trailing singleton dimension).
@inline get_step(var::AbstractArray{T, 1}, step::Integer) where {T} = view(var, :, step)
@inline get_step(var::AbstractArray{T, 2}, step::Integer) where {T} = view(var, :, step)
@inline get_step(var::AbstractArray{T, 3}, step::Integer) where {T} = view(var, :, :, step)

# LOWER TRIANGULAR ARRAYS
# `step` always indexes the time dimension `T` in the variable's ArrayDimensions (LMT, LMZT).
# Variables without a time dimension (LM, LMZ) have no step dimension to index into, so the
# full array is returned as a view and `step` is ignored (it should be 1).

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the time dimension `T` in the variable's `ArrayDimensions`, for prognostic
variables e.g. used for the leapfrog time step. These methods are for spectral variables
WITHOUT a time dimension (`LM`, `LMZ`), which have no step dimension to select from: the full
array is returned as a view and `step` is ignored (it should be 1)."""
@inline get_step(var::LowerTriangularArray{T, 1, A, S, ArrayDimensions.LM}, step::Integer) where {T, A, S} = lta_view(var, :)
@inline get_step(var::LowerTriangularArray{T, 2, A, S, ArrayDimensions.LMZ}, step::Integer) where {T, A, S} = lta_view(var, :, :)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the time dimension `T` in the variable's `ArrayDimensions`, for prognostic
variables e.g. used for the leapfrog time step. These methods are for spectral variables WITH a
time dimension (`LMT` horizontal + time, `LMZT` horizontal + vertical + time), whose last
dimension is the step dimension."""
@inline get_step(var::LowerTriangularArray{T, 2, A, S, ArrayDimensions.LMT}, step::Integer) where {T, A, S} = lta_view(var, :, step)
@inline get_step(var::LowerTriangularArray{T, 3, A, S, ArrayDimensions.LMZT}, step::Integer) where {T, A, S} = lta_view(var, :, :, step)

# FIELDS
# same as for LowerTriangularArrays: `step` always indexes the time dimension `T` in the
# variable's ArrayDimensions (XYT, XYZT), variables without one (XY, XYZ) return the full view.

"""$(TYPEDSIGNATURES)
Get the i-th step of a field as a view (wrapped into the same type as the input variable).
"step" refers to the time dimension `T` in the variable's `ArrayDimensions`, for prognostic
variables e.g. used for the leapfrog time step. These methods are for fields WITHOUT a time
dimension (`XY`, `XYZ`), which have no step dimension to select from: the full array is
returned as a view and `step` is ignored (it should be 1)."""
@inline get_step(var::AbstractField{T, 1, A, G, ArrayDimensions.XY}, step::Integer) where {T, A, G} = field_view(var, :)
@inline get_step(var::AbstractField{T, 2, A, G, ArrayDimensions.XYZ}, step::Integer) where {T, A, G} = field_view(var, :, :)

"""$(TYPEDSIGNATURES)
Get the i-th step of a field as a view (wrapped into the same type as the input variable).
"step" refers to the time dimension `T` in the variable's `ArrayDimensions`, for prognostic
variables e.g. used for the leapfrog time step. These methods are for fields WITH a time
dimension (`XYT` horizontal + time, `XYZT` horizontal + vertical + time), whose last dimension
is the step dimension."""
@inline get_step(var::AbstractField{T, 1, A, G, ArrayDimensions.XYT}, step::Integer) where {T, A, G} = field_view(var, :, step)
@inline get_step(var::AbstractField{T, 2, A, G, ArrayDimensions.XYT}, step::Integer) where {T, A, G} = field_view(var, :, step)
@inline get_step(var::AbstractField{T, 3, A, G, ArrayDimensions.XYZT}, step::Integer) where {T, A, G} = field_view(var, :, :, step)

# anything that can decide which variable step to get
const STEP_COMPONENT = Union{AbstractModelComponent, SpeedyTransforms.AbstractSpectralTransform}

"""Time steppers control on which step (index of the step dimension, typically used to store two time steps, e.g. t-dt, t)
a variable is evaluated for a given term or in a given computation of the dynamical core.
Two types of variables have step dimensions: Prognostic variables and tendencies.
Leapfrog uses 2 steps for the prognostic variables but one tendency;
Adams-Bashforth is a multi-step method storing one step for the prognostic variables and multiple steps for the tendencies.

This is implemented in the dynamical core via calling `get_step`, and particularly
`get_prognostic_step`, `get_tendency_step` which dispatch on the time stepper and the component for which the step is needed.
Dispatch in these functions has to be

    var::Any, ::AbstractTimeStepper, ::STEP_COMPONENT, ::AbstractModel

whereby the model can be left out. Default is to return step index 1.
This means that every time stepper can define which step to call in the `get_..._step` methods
by implementing `which_step`, `which_prognostic_step`, or `which_tendency_step`.
With that function signature. So dispatch allows to distinguish the step between model components,
between models but also whether a variable is spectral or a grid variable."""
which_step

# methods independent of model
@inline get_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = get_step(var, which_step(var, TS, C))
@inline get_prognostic_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = get_step(var, which_prognostic_step(var, TS, C))
@inline get_tendency_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = get_step(var, which_tendency_step(var, TS, C))

# methods dependent on model
@inline get_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, M::AbstractModel) = get_step(var, which_step(var, TS, C, M))
@inline get_prognostic_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, M::AbstractModel) = get_step(var, which_prognostic_step(var, TS, C, M))
@inline get_tendency_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, M::AbstractModel) = get_step(var, which_tendency_step(var, TS, C, M))

# if dispatch over model is not defined then fallback to dispatch without model
@inline which_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, ::AbstractModel) = which_step(var, TS, C)
@inline which_prognostic_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, ::AbstractModel) = which_prognostic_step(var, TS, C)
@inline which_tendency_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT, ::AbstractModel) = which_tendency_step(var, TS, C)

# fallback to 1 if not defined
@inline which_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = 1
@inline which_prognostic_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = 1
@inline which_tendency_step(var, TS::AbstractTimeStepper, C::STEP_COMPONENT) = 1
