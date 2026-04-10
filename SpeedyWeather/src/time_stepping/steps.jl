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
as defined by `get_step`. "Steps" refer to the last dimension, for prognostic variables e.g. used for the leapfrog time step."""
get_steps

get_steps(var::AbstractArray{T, 1}) where {T} = (var,)
get_steps(var::AbstractArray{T, 2}) where {T} = ntuple(step -> get_step(var, step), size(var, 2))
get_steps(var::AbstractArray{T, 3}) where {T} = ntuple(step -> get_step(var, step), size(var, 3))

export get_step

# LowerTriangularArrays
# for 2D spectral variables step can be 1 that'll be the ignored additional singleton dimension
# otherwise an error is thrown
@inline get_step(var::LowerTriangularArray{T, 1}, step::Integer) where {T} = lta_view(var, :, step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 2D spectral variable (horizontal only) with steps in the 3rd dimension."""
@inline get_step(var::LowerTriangularArray{T, 2}, step::Integer) where {T} = lta_view(var, :, step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a LowerTriangularArray as a view (wrapped into a LowerTriangularArray).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 3D spectral variable (horizontal + vertical) with steps in the 4rd dimension."""
@inline get_step(var::LowerTriangularArray{T, 3}, step::Integer) where {T} = lta_view(var, :, :, step)

# FIELDS
# for 2D fields step can be 1 that'll be the ignored additional singleton dimension
# otherwise an error is thrown
@inline get_step(var::AbstractField{T, 1}, step::Integer) where {T} = field_view(var, :, step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a 3D field as a view (wrapped into the same type as the input variable).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 2D field (horizontal only) with steps in the 3rd dimension."""
@inline get_step(var::AbstractField{T, 2}, step::Integer) where {T} = field_view(var, :, step)

"""$(TYPEDSIGNATURES)
Get the i-th step of a 4D field as a view (wrapped into the same type as the input variable).
"step" refers to the last dimension, for prognostic variables e.g. used for the leapfrog time step.
This method is for a 3D field (horizontal + vertical) with steps in the 4rd dimension."""
@inline get_step(var::AbstractField{T, 3}, step::Integer) where {T} = field_view(var, :, :, step)

# get_step depending on time stepping method
@inline get_step(var, ::AbstractTimeStepper) = get_step(var, 1)
@inline get_step(var, ::AbstractTimeStepper, ::AbstractModelComponent) = get_step(var, 1)

@inline get_prognostic_step(var, TS::AbstractTimeStepper, args...) = get_step(var, TS, args...)
@inline get_tendency_step(var, TS::AbstractTimeStepper, args...) = get_step(var, TS, args...)