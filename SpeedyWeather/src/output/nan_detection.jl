abstract type AbstractNaNDetector <: AbstractModelComponent end

export NaNDetector

"""
Detects NaNs (or Infs) in a simulation, which stops it. NaNs propagate through the spectral transform
to all coefficients, so whether the 1-0 mode of the surface vorticity is finite is recorded on every
time step (without synchronizing a GPU), together with the time step and the simulation time. The record
is checked on the CPU every `interval` time steps and at the end, the warning then reports the time step
and simulation time at which NaNs first occurred, also stored in `step` and `time`. Used as
`Feedback(nan_detection = NaNDetector())`, `nan_detection = nothing` disables it.
$(TYPEDFIELDS)"""
@kwdef mutable struct NaNDetector <: AbstractNaNDetector
    "[OPTION] Interval in time steps between checks of the record (copying it to the CPU)"
    interval::Int = 50

    "[DERIVED] were NaNs detected?"
    detected::Bool = false

    "[DERIVED] time step at which NaNs were first detected, 0 if none"
    step::Int = 0

    "[DERIVED] simulation time at which NaNs were first detected"
    time::DateTime = DateTime(2000)

    "[DERIVED] time steps recorded since the last check"
    recorded_steps::Vector{Int} = zeros(Int, interval)

    "[DERIVED] simulation times of the time steps recorded since the last check"
    recorded_times::Vector{DateTime} = fill(DateTime(2000), interval)

    "[DERIVED] number of time steps recorded since the last check"
    nrecorded::Int = 0
end

# the record on the device, 1 if the 1-0 mode of the surface vorticity was not finite in a recorded time step
variables(detector::NaNDetector) = (
    DynamicsVariable(
        :nan_record, VectorDim(max(1, detector.interval)), namespace = :nan_detection,
        desc = "1 if the 1-0 mode of the surface vorticity was NaN/Inf in a recorded time step, 0 otherwise", units = "1"
    ),
)

function initialize!(detector::NaNDetector, vars::Variables, model::AbstractModel)
    interval = length(vars.dynamics.nan_detection.nan_record)
    detector.detected = false
    detector.step = 0
    detector.recorded_steps = zeros(Int, interval)
    detector.recorded_times = fill(vars.prognostic.clock.time, interval)
    detector.nrecorded = 0
    return nothing
end

"""$(TYPEDSIGNATURES)
Were NaNs detected? `false` without NaN detection (`nothing`)."""
nans_detected(detector::NaNDetector) = detector.detected
nans_detected(::Nothing) = false

"""$(TYPEDSIGNATURES)
Record whether the 1-0 mode of the surface vorticity is finite at time step `step`, and check the
record on the CPU every `interval` recorded time steps or if `check = true` (e.g. on the last time step).
Warns with the time step and simulation time at which NaNs first occurred."""
function detect_nans!(detector::NaNDetector, vars::Variables, model::AbstractModel, step::Integer; check::Bool = false)
    detector.detected && return nothing         # nothing to do once NaNs were detected
    record = vars.dynamics.nan_detection.nan_record
    k = detector.nrecorded += 1
    vor = get_prognostic_step(vars.prognostic.vorticity, model.time_stepping, detector)
    record_nan!(record, k, vor.data, size(vor.data, 2), model.architecture)
    detector.recorded_steps[k] = step
    detector.recorded_times[k] = vars.prognostic.clock.time

    (check || k == length(record)) || return nothing

    # check the record on the CPU (synchronizes a GPU)
    detector.nrecorded = 0
    first_nan = findfirst(!iszero, Array(view(record, 1:k)))
    isnothing(first_nan) && return nothing
    detector.detected = true
    detector.step = detector.recorded_steps[first_nan]
    detector.time = detector.recorded_times[first_nan]
    @warn "NaN or Inf detected at time step $(detector.step) ($(detector.time))"
    return nothing
end

detect_nans!(::Nothing, args...; kwargs...) = nothing

# record[k] = !isfinite(vor[2, layer]) on the device, without synchronizing a GPU
function record_nan!(record::AbstractVector, k::Integer, vor::AbstractArray, layer::Integer, ::AbstractArchitecture)
    view(record, k:k) .= .!isfinite.(view(vor, 2, layer:layer))
    return nothing
end

# Reactant arrays outside a compiled program (where the feedback happens): scalar indexing instead
function record_nan!(record::AbstractVector, k::Integer, vor::AbstractArray, layer::Integer, ::ReactantDevice)
    GPUArrays.@allowscalar record[k] = !isfinite(vor[2, layer])
    return nothing
end
