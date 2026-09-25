export EulerForward

"""Euler forward time stepping, used by default for the ocean and land namespaces
with `Leapfrog` (see its `ocean` and `land` fields). Defined by the following fields
$(TYPEDFIELDS)"""
mutable struct EulerForward{NF, S, B, MS} <: AbstractTimeStepper
    "[OPTION] Time step for T32, scale linearly to spectral resolution `truncation`"
    Δt_at_T32::S

    "[OPTION] Adjust `Δt_at_T32` with the output `interval` to output exactly after integer time steps"
    adjust_with_output::B

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt::NF
end

"""$(TYPEDSIGNATURES)
Generator function for an EulerForward struct using `spectral_grid`
for the resolution information."""
function EulerForward(
        spectral_grid::SpectralGrid;
        Δt_at_T32 = Minute(40),
        adjust_with_output = true,
    )
    (; NF, truncation) = spectral_grid
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T32), truncation, DEFAULT_RADIUS, adjust_with_output)
    Δt::NF = Δt_millisec.value / 1000
    return EulerForward(Second(Δt_at_T32), adjust_with_output, Δt_millisec, Δt)
end

"""$(TYPEDSIGNATURES) Initialize `E` by recalculating the time step for the resolution of `model`."""
function initialize!(E::EulerForward, model::AbstractModel)
    calculate_Δt!(E, model)
    return nothing
end

# 1 prognostic and 1 tendency step, and read step 1, are the fallbacks in steps.jl

"""$(TYPEDSIGNATURES)
Euler forward step `var += Δt/scale * tendency`."""
function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        clock::Clock,
        time_stepping::EulerForward,
        implicit,
        ::AbstractModel,
        scale::Real = 1,
    )
    Δt = time_step(time_stepping, clock)
    Δt /= oftype(Δt, scale)
    @boundscheck size(var) == size(tendency) || throw(BoundsError())
    launch!(architecture(var), LinearWorkOrder, size(var), euler_forward_kernel!, var, tendency, Δt)
    return nothing
end

@kernel inbounds = true function euler_forward_kernel!(var, tendency, Δt)
    ij = @index(Global, Linear)
    var[ij] += Δt * tendency[ij]
end
