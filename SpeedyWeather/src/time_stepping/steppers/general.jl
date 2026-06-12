"""($TYPEDSIGNATURES) Extend in case a time stepper requires spin up steps,
e.g. Leapfrog starts with 1 Euler step that isn't counted for the clock but requires
one more step in the main time loop."""
spin_up_steps(::AbstractTimeStepper) = 0
default_time_step(L::AbstractTimeStepper) = L.Δt

"""$(TYPEDSIGNATURES)
Computes the time step in [ms]. `Δt_at_T31` is always scaled with the resolution `trunc` 
of the model. In case `adjust_Δt_with_output` is true, the `Δt_at_T31` is additionally 
adjusted to the closest divisor of `interval` so that the output time axis is keeping
`interval` exactly."""
function get_Δt_millisec(
        Δt_at_T31::TimePeriod,
        trunc,
        radius,
        adjust_with_output::Bool,
        interval::TimePeriod = DEFAULT_OUTPUT_INTERVAL,
    )
    # linearly scale Δt with trunc+1 (which are often powers of two)
    resolution_factor = (DEFAULT_TRUNC + 1) / (trunc + 1)

    # radius also affects grid spacing, scale proportionally
    radius_factor = radius / DEFAULT_RADIUS

    # maybe rename to _at_trunc_and_radius?
    Δt_at_trunc = Second(Δt_at_T31).value * resolution_factor * radius_factor

    if adjust_with_output && (interval > Millisecond(0))
        k = round(Int, Second(interval).value / Δt_at_trunc)
        divisors = Primes.divisors(Millisecond(interval).value)
        sort!(divisors)
        i = findfirst(x -> x >= k, divisors)
        k_new = isnothing(i) ? k : divisors[i]
        Δt_millisec = Millisecond(round(Int, Millisecond(interval).value / k_new))

        # provide info when time step is significantly shortened or lengthened
        Δt_millisec_unadjusted = round(Int, 1000 * Δt_at_trunc)
        Δt_ratio = Δt_millisec.value / Δt_millisec_unadjusted

        if abs(Δt_ratio - 1) > 0.05     # print info only when +-5% changes
            p = round(Int, (Δt_ratio - 1) * 100)
            ps = p > 0 ? "+" : ""
            @info "Time step changed from $Δt_millisec_unadjusted to $Δt_millisec ($ps$p%) to match output frequency."
        end
    else
        Δt_millisec = Millisecond(round(Int, 1000 * Δt_at_trunc))
    end

    return Δt_millisec
end

"""$(TYPEDSIGNATURES)
Change time step of timestepper `L` to `Δt` (unscaled)
and disables adjustment to output frequency.
`set!` can be used before or after `initialize!(model)`."""
function set!(
        L::AbstractTimeStepper,
        Δt::Period;                         # unscaled time step in Second, Minute, ...
        radius = L.Δt_sec/L.Δt              # get radius from scaled time step
    )
    # if set! is used before `initialize!(model)` then the recalculation of
    # Δt_at_T31 will make sure the desired time step is used when initialize!(model.time_stepping, ...) happens
    # if set! is used after `initialize!(model)` then all fields are set consistently

    # get truncation/resolution factor implicitly from Δt at this resolution vs default T31 resolution
    resolution_factor = L.Δt_sec / Second(L.Δt_at_T31).value

    L.Δt_millisec = Millisecond(Δt)         # recalculate all Δt fields
    L.Δt_sec = Millisecond(L.Δt_millisec).value / 1000
    L.Δt = L.Δt_sec / radius

    # recalculate the default time step at resolution T31 to be consistent
    L.Δt_at_T31 = Second(round(Int, L.Δt_sec * resolution_factor))

    # given Δt was manually set disallow adjustment to output frequency
    L.adjust_with_output = false
    return L
end

# also allow for keyword arguments
set!(L::AbstractTimeStepper; Δt::Period, kwargs...) = set!(L, Δt; kwargs...)

function calculate_Δt!(L::AbstractTimeStepper, model::AbstractModel)
    (; trunc) = model.spectral_grid
    (; radius) = model.planet
    interval = get_interval(model.output)

    # take radius from planet and recalculate time step and possibly adjust with output dt
    L.Δt_millisec = get_Δt_millisec(L.Δt_at_T31, trunc, radius, L.adjust_with_output, interval)
    L.Δt_sec = L.Δt_millisec.value / 1000
    L.Δt = L.Δt_sec / radius

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
