"""
$(TYPEDSIGNATURES)
Linear virtual temperature for `model::PrimitiveDry`: Just copy over
arrays from `temp` to `temp_virt` at timestep `lf` in spectral space
as humidity is zero in this `model`."""
function linear_virtual_temperature!(
        vars::Variables,
        model::PrimitiveDry,
    )
    Tᵥ = vars.dynamics.virtual_temperature

    # For Leapfrog this term has to be evaluted on the previous time step
    # as the implicit corrections will move it to the current as done for 
    # all linear gravity-wave related terms, just denote this with `LinearDynamicalCore`
    # here, the time stepper then decides which step to return
    T = get_prognostic_step(vars.prognostic.temperature, model.time_stepping, LinearDynamicalCore())
    Tᵥ .= T
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculates a linearised virtual temperature Tᵥ as

    Tᵥ = T + Tₖμq

With absolute temperature T, layer-average temperarture Tₖ (computed in `temperature_average!`),
specific humidity q and

    μ = (1-ξ)/ξ, ξ = R_dry/R_vapor.

in spectral space."""
function linear_virtual_temperature!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    Tᵥ = vars.dynamics.virtual_temperature
    μ = model.atmosphere.μ_virt_temp
    Tₖ = vars.dynamics.average_temperature_profile

    # For Leapfrog this term has to be evaluted on the previous time step
    # as the implicit corrections will move it to the current as done for 
    # all linear gravity-wave related terms, just denote this with `LinearDynamicalCore`
    # here, the time stepper then decides which step to return
    T = get_prognostic_step(vars.prognostic.temperature, model.time_stepping, LinearDynamicalCore())
    q = get_prognostic_step(vars.prognostic.humidity, model.time_stepping, LinearDynamicalCore())

    # TODO check that doing a non-linear virtual temperature in grid-point space
    # but a linear virtual temperature in spectral space to avoid another transform
    # does not cause any problems. Alternative do the transform or have a linear
    # virtual temperature in both grid and spectral space

    # TODO: broadcast with LTA doesn't work here becasue of a broadcast conflict
    # (Tₖ and q are different dimensions and array types)
    @. Tᵥ.data = T.data + (Tₖ' * μ) * q.data
    return nothing
end

@inline virtual_temperature(T, q, A::AbstractWetAtmosphere) = virtual_temperature(T, q, A.μ_virt_temp)
@inline linear_virtual_temperature(T, q, Tₖ, A::AbstractWetAtmosphere) = linear_virtual_temperature(T, q, Tₖ, A.μ_virt_temp)
@inline absolute_temperature(Tᵥ, q, A::AbstractWetAtmosphere) = absolute_temperature(Tᵥ, q, A.μ_virt_temp)

@inline virtual_temperature(T, q, μ) = T * (1 + μ * q)
@inline linear_virtual_temperature(T, q, Tₖ, μ) = T + (Tₖ * μ) * q
@inline absolute_temperature(Tᵥ, q, μ) = Tᵥ / (1 + μ * q)

# For dry atmospheres, virtual temperature is just (absolute) temperature
@inline virtual_temperature(T, q, ::AbstractDryAtmosphere) = T
@inline linear_virtual_temperature(T, q, Tₖ, ::AbstractDryAtmosphere) = T
@inline absolute_temperature(T, q, ::AbstractDryAtmosphere) = T
