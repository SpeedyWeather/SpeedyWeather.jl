"""($TYPEDSIGNATURES) Extend in case a time stepper requires spin up steps,
e.g. Leapfrog starts with 1 Euler step that isn't counted for the clock but requires
one more step in the main time loop."""
spin_up_steps(::AbstractTimeStepper) = 0
default_time_step(L::AbstractTimeStepper) = L.Δt

"""$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T32` is always scaled with the resolution `truncation` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T32` is additionally 
adjusted to the closest divisor of `interval` so that the output time axis is keeping
`interval` exactly."""
function get_Δt_millisec(
        Δt_at_T32::TimePeriod,
        truncation,                     # spectral truncation
        radius,                         # of the planet (determines Δx)
        adjust_with_output::Bool,
        interval::TimePeriod = DEFAULT_OUTPUT_INTERVAL,
    )
    resolution_factor = DEFAULT_TRUNCATION / truncation     # linearly scale Δt with truncation
    radius_factor = radius / DEFAULT_RADIUS                 # radius also affects grid spacing, scale proportionally
    Δt_scaled = Second(Δt_at_T32).value * resolution_factor * radius_factor

    if adjust_with_output && (interval > Millisecond(0))
        k = round(Int, Second(interval).value / Δt_scaled)
        divisors = Primes.divisors(Millisecond(interval).value)
        sort!(divisors)
        i = findfirst(x -> x >= k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(interval).value / k_new))

        # provide info when time step is significantly shortened or lengthened
        Δt_millisec_unadjusted = round(Int, 1000 * Δt_scaled)
        Δt_ratio = Δt_millisec.value / Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05     # print info only when +-5% changes
            p = round(Int, (Δt_ratio - 1) * 100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else
        Δt_millisec = Millisecond(round(Int, 1000 * Δt_scaled))
    end

    return Δt_millisec
end

"""$(TYPEDSIGNATURES)
Factor between the time step at T32 and the time step at resolution `trunc` on a planet of
`radius`, i.e. `Δt = Δt_at_T32 * resolution_factor`, matching [`get_Δt_millisec`](@ref)."""
resolution_factor(truncation::Integer, radius::Real) =
    DEFAULT_TRUNCATION / truncation * radius / DEFAULT_RADIUS

"""$(TYPEDSIGNATURES) Resolution factor of `model`, from its `trunc` and planetary radius."""
resolution_factor(model::AbstractModel) =
    resolution_factor(model.spectral_grid.truncation, model.planet.radius)

"""$(TYPEDSIGNATURES)
Resolution factor inferred from a time stepper's own state, `Δt/Δt_at_T32`. Approximate: if the
current `Δt` was adjusted to the output frequency (`adjust_with_output`) then that adjustment is
folded into the ratio. Pass the `model` to [`set!`](@ref) to use the exact factor instead."""
resolution_factor(L::AbstractTimeStepper) = L.Δt / Second(L.Δt_at_T32).value

"""$(TYPEDSIGNATURES)
Change the time step of timestepper `L` to `Δt` (unscaled) and disable adjustment to the output
frequency. `set!` can be used before or after `initialize!(model)`: `Δt_at_T32` is the quantity
`calculate_Δt!` derives the time step from, so it is back-computed here with `resolution_factor`
to keep a later re-initialization at the requested `Δt`.

Without a `model` the resolution factor is inferred from `L` itself and is only approximate when
the current time step had been adjusted to the output frequency; pass `model` (or use
`set!(model, Δt=...)`) for an exact round trip."""
function set!(
        L::AbstractTimeStepper,
        Δt::Period,
        factor::Real = resolution_factor(L),
    )
    L.Δt_millisec = Millisecond(Δt)         # recalculate all Δt fields
    L.Δt = Millisecond(L.Δt_millisec).value / 1000

    # Δt = Δt_at_T32 * factor (see get_Δt_millisec), so invert to keep Δt_at_T32 consistent.
    # Stored in whole seconds, so a re-initialization can be off by sub-second rounding.
    L.Δt_at_T32 = Second(round(Int, L.Δt / factor))

    # given Δt was manually set disallow adjustment to output frequency
    L.adjust_with_output = false
    return L
end

"""$(TYPEDSIGNATURES) Set the time step of `L` using the exact resolution factor of `model`."""
set!(L::AbstractTimeStepper, Δt::Period, model::AbstractModel) =
    set!(L, Δt, resolution_factor(model))

# also allow for keyword arguments
set!(L::AbstractTimeStepper; Δt::Period) = set!(L, Δt)
set!(L::AbstractTimeStepper, model::AbstractModel; Δt::Period) = set!(L, Δt, model)

function calculate_Δt!(L::AbstractTimeStepper, model::AbstractModel)
    (; truncation) = model.spectral_grid
    (; radius) = model.planet
    interval = get_interval(model.output)

    # take radius from planet and recalculate time step and possibly adjust with output dt
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T32, truncation, radius, L.adjust_with_output, interval)
    L.Δt = L.Δt_millisec.value / 1000

    # check how time steps from time integration and output align
    if L.adjust_with_output
        n = round(Int, Millisecond(interval).value / L.Δt_millisec.value)
        nΔt = n * L.Δt_millisec
        if nΔt != interval
            @warn "$n steps of Δt = $(L.Δt_millisec.value)ms yield output every " *
                "$(nΔt.value)ms (=$(nΔt.value / 1000)s), but interval = $(interval.value)s"
        end
    end
    return nothing
end

# extend in case time steppers want to initialize variables differently before the first
# simulation is run, e.g. Leapfrog may want to copy prognostic steps
@inline initialize!(::Variables, ::AbstractTimeStepper, ::AbstractModel) = nothing
