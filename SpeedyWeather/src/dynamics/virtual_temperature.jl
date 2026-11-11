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
