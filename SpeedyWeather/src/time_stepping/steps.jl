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
as defined by `get_step`. "Steps" refer to the last dimension, for prognostic variables
e.g. used for the leapfrog time step."""
get_steps

get_steps(var::AbstractArray{T, 1}) where {T} = (var,)
get_steps(var::AbstractArray{T, 2}) where {T} = ntuple(step -> get_step(var, step), size(var, 2))
get_steps(var::AbstractArray{T, 3}) where {T} = ntuple(step -> get_step(var, step), size(var, 3))

export get_step

"""$(TYPEDSIGNATURES)
Select step dimension from variable, when no step as 2nd argument provided select las index
as this typically presents the "current" step (and not any previous ones). But this depends
on the time stepping a variable with step dimension was created for."""
get_step(var) = get_step(var, size(var, ndims(var)))

# Plain Arrays
# Inside GPU kernels `Adapt.adapt_structure(to, field::AbstractField) = adapt(to, field.data)`
# unwraps a Field to its bare device array, so `get_step` must also work on plain arrays
# (otherwise the device-side MethodError aborts kernel compilation). Same semantics as the
# Field/LowerTriangularArray methods below: `step` indexes the last dimension; for arrays
# without an explicit step dimension `step` must be 1 (trailing singleton dimension).
@inline get_step(var::AbstractArray{T, 1}, step::Integer) where {T} = view(var, :, step)
@inline get_step(var::AbstractArray{T, 2}, step::Integer) where {T} = view(var, :, step)
@inline get_step(var::AbstractArray{T, 3}, step::Integer) where {T} = view(var, :, :, step)

# A view built from a *runtime* last-dimension offset cannot be type-analysed by Enzyme on Julia
# ≥ 1.11 once that view participates in a differentiated KernelAbstractions kernel — it throws an
# `EnzymeNoTypeError`. A view with a *compile-time* offset is handled natively and correctly (the
# aliasing read/write semantics that the leapfrog relies on are preserved). `_literal_step` therefore
# resolves the small, statically-bounded step index to a literal via a short branch ladder, so every
# step view carries a compile-time offset.
#
# Note: a custom `get_step` `EnzymeRules` rule was tried instead (to keep this branchless for Reactant)
# but rejected — `get_step` returns a view used *mutably* (read and written, aliasing the prognostic)
# in the leapfrog, which a custom rule cannot model; `EnzymeTestUtils.test_reverse` confirmed it gave
# wrong gradients. The branch is resolved at trace time whenever `step` is a concrete value (the usual
# case, since the clock step counter is a scalar), so it does not introduce control flow under Reactant.
@inline function _literal_step(builder::F, step::Integer) where {F}
    step == 1 && return builder(1)
    step == 2 && return builder(2)
    return builder(step)
end

# LowerTriangularArrays
# for 2D spectral variables step can be 1 that'll be the ignored additional singleton dimension
# otherwise an error is thrown
@inline get_step(var::LowerTriangularArray{T, 1}, step::Integer) where {T} = _literal_step(d -> lta_view(var, :, d), step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 2D spectral variable (horizontal only) with steps in the 3rd dimension."""
@inline get_step(var::LowerTriangularArray{T, 2}, step::Integer) where {T} = _literal_step(d -> lta_view(var, :, d), step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 3D spectral variable (horizontal + vertical) with steps in the 4rd dimension."""
@inline get_step(var::LowerTriangularArray{T, 3}, step::Integer) where {T} = _literal_step(d -> lta_view(var, :, :, d), step)

# FIELDS
# for 2D fields step can be 1 that'll be the ignored additional singleton dimension
# otherwise an error is thrown
@inline get_step(var::AbstractField{T, 1}, step::Integer) where {T} = _literal_step(d -> field_view(var, :, d), step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a 3D field as a view (wrapped into the same type as the input variable).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 2D field (horizontal only) with steps in the 3rd dimension."""
@inline get_step(var::AbstractField{T, 2}, step::Integer) where {T} = _literal_step(d -> field_view(var, :, d), step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a 4D field as a view (wrapped into the same type as the input variable).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 3D field (horizontal + vertical) with steps in the 4rd dimension."""
@inline get_step(var::AbstractField{T, 3}, step::Integer) where {T} = _literal_step(d -> field_view(var, :, :, d), step)

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
