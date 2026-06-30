abstract type AbstractLongwaveTransmissivity <: AbstractLongwave end

export TransparentLongwaveTransmissivity
TransparentLongwaveTransmissivity(SG::SpectralGrid) = ConstantLongwaveTransmissivity(SG, transmissivity = 1)

export ConstantLongwaveTransmissivity
@parameterized @kwdef struct ConstantLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    @param transmissivity::NF = 0.6 (bounds = 0 .. 1,)
end
Adapt.@adapt_structure ConstantLongwaveTransmissivity
ConstantLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = ConstantLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::ConstantLongwaveTransmissivity, ::AbstractModel) = nothing
@propagate_inbounds function transmissivity!(ij, vars, CLT::ConstantLongwaveTransmissivity, model)
    t = vars.scratch.grid.a
    nlayers = size(t, 2)

    τ = -log(CLT.transmissivity)            # total optical depth of the atmosphere
    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]

    for k in 1:nlayers
        Δσₖ = pressure_thickness(k, pₛ, coord) / pₛ   
        t[ij, k] = exp(-τ * Δσₖ)            # transmissivity through layer k
    end
    return t
end

export FriersonLongwaveTransmissivity
@parameterized @kwdef struct FriersonLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    "[OPTION] Optical depth at the equator"
    @param τ₀_equator::NF = 6 (bounds = Nonnegative,)

    "[OPTION] Optical depth at the poles"
    @param τ₀_pole::NF = 1.5 (bounds = Nonnegative,)

    "[OPTION] Fraction to mix linear and quadratic profile"
    @param fₗ::NF = 0.1 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure FriersonLongwaveTransmissivity
FriersonLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = FriersonLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::FriersonLongwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(ij, vars, transmissivity::FriersonLongwaveTransmissivity, model)

    # use scratch array to compute transmissivity t
    t = vars.scratch.grid.a
    nlayers = size(t, 2)
    NF = eltype(t)

    # but the longwave optical depth follows some idealised humidity lat-vert distribution
    (; τ₀_equator, τ₀_pole, fₗ) = transmissivity

    # coordinates
    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]
    j = model.geometry.whichring[ij]
    sinlat = model.geometry.sinlat[j]   # sin(latitude), avoids calling sind() on traced scalars

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    τ_above::NF = 0
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2
    for k in 1:nlayers              # loop over half levels below
        σₖ = pressure_below(k, pₛ, coord) / pₛ
        τ_below = τ₀ * (fₗ * σₖ + (1 - fₗ) * σₖ^4)
        t[ij, k] = exp(-(τ_below - τ_above))
        τ_above = τ_below
    end

    # return so the radiative_trasfer uses the right scratch array
    return t
end

"""$(TYPEDSIGNATURES)
Precompute the (state-independent) longwave transmissivity once into
a `nlayers × nlat` matrix `t` (one column per latitude ring). Computing on the
host with libm `exp` and transferring the SAME array to the device makes the
transmissivity bit-identical across CPU/GPU/Reactant."""
function fill_longwave_transmissivity!(t::AbstractMatrix, CLT::ConstantLongwaveTransmissivity, model)
    nlayers, nlat = size(t)
    # pull geometry to host (Vector) so the loop runs host libm `exp`, not device/XLA
    # `exp`; on CPU this is a no-op, under Reactant it transfers ConcretePJRTArray → Vector
    dσ = on_architecture(CPU(), model.geometry.σ_levels_thick)
    τ = -log(CLT.transmissivity)            # total optical depth of the atmosphere
    for j in 1:nlat, k in 1:nlayers
        t[k, j] = exp(-τ * dσ[k])           # transmissivity through layer k (lat-independent)
    end
    return t
end

function fill_longwave_transmissivity!(t::AbstractMatrix, transmissivity::FriersonLongwaveTransmissivity, model)
    nlayers, nlat = size(t)
    NF = eltype(t)
    (; τ₀_equator, τ₀_pole, fₗ) = transmissivity
    # pull geometry to host (Vector) so the loop runs host libm `exp`, not device/XLA `exp`
    σ = on_architecture(CPU(), model.geometry.σ_levels_half)
    sinlat = on_architecture(CPU(), model.geometry.sinlat)          # sin(latitude) per ring
    for j in 1:nlat
        τ_above = zero(NF)
        τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat[j]^2
        for k in 2:(nlayers + 1)
            τ_below = τ₀ * (fₗ * σ[k] + (1 - fₗ) * σ[k]^4)
            t[k - 1, j] = exp(-(τ_below - τ_above))
            τ_above = τ_below
        end
    end
    return t
end
