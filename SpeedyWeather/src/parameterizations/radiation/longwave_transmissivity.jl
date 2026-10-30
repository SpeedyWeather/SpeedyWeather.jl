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
    θ = model.geometry.latds[ij]

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    τ_above::NF = 0

    # TODO: Replace `sin(deg2rad(θ))` with `sind(θ)` once JuliaGPU/AMDGPU.jl#1041 
    # is merged/released and `sind` is supported on AMD GPUs.k
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sin(deg2rad(θ))^2
    for k in 1:nlayers              # loop over half levels below
        σₖ = pressure_below(k, pₛ, coord) / pₛ
        τ_below = τ₀ * (fₗ * σₖ + (1 - fₗ) * σₖ^4)
        t[ij, k] = exp(-(τ_below - τ_above))
        τ_above = τ_below
    end

    # return so the radiative_trasfer uses the right scratch array
    return t
end

export ByrneOGormanLongwaveTransmissivity
"""Longwave transmissivity with an optical depth that depends on specific humidity and CO₂, following
Byrne and O'Gorman, 2013, J. Climate, https://doi.org/10.1175/JCLI-D-12-00262.1, as also used in Isca.
The optical depth of layer k is `dτ = (a + c*log(CO₂/CO₂ᵣ) + b*q) * Δp / p₀` with pressure thickness `Δp`,
specific humidity `q` [kg/kg], reference pressure `p₀` from `model.atmosphere`, a dry absorption `a` (well-mixed greenhouse gases
at the reference concentration `CO₂ᵣ`), a CO₂ sensitivity `c` and a water vapor absorption `b`.
CO₂ [ppm] is read from `model.greenhouse_gases.co2` if present, otherwise `CO₂ᵣ` is used.
In contrast to the `FriersonLongwaveTransmissivity` there is no prescribed latitude-vertical profile,
the optical depth follows the model's humidity, introducing a water vapor feedback. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ByrneOGormanLongwaveTransmissivity{NF} <: AbstractLongwaveTransmissivity
    "[OPTION] Optical depth of the dry atmosphere (well-mixed greenhouse gases) per p₀ at `reference_co2` [1]"
    @param dry_absorption::NF = 0.8678 (bounds = Nonnegative,)

    "[OPTION] Optical depth from water vapor per p₀ and per specific humidity [1/(kg/kg)]"
    @param water_vapor_absorption::NF = 1997.9 (bounds = Nonnegative,)

    "[OPTION] Optical depth per p₀ per e-fold of CO₂ relative to `reference_co2` [1]"
    @param co2_sensitivity::NF = 0.2023 (bounds = Nonnegative,)

    "[OPTION] CO₂ concentration at which the dry optical depth is `dry_absorption` [ppm]"
    reference_co2::NF = 360
end

Adapt.@adapt_structure ByrneOGormanLongwaveTransmissivity
ByrneOGormanLongwaveTransmissivity(SG::SpectralGrid; kwargs...) = ByrneOGormanLongwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::ByrneOGormanLongwaveTransmissivity, ::AbstractModel) = nothing

# CO₂ [ppm] from the greenhouse gases if present, otherwise the reference
co2_concentration(prognostic, reference_co2) =
    haskey(prognostic, :greenhouse_gases) && haskey(prognostic.greenhouse_gases, :co2) ?
    convert(typeof(reference_co2), prognostic.greenhouse_gases.co2[]) : reference_co2

@propagate_inbounds function transmissivity!(ij, vars, transmissivity::ByrneOGormanLongwaveTransmissivity, model)
    t = vars.scratch.grid.a                         # use scratch array to compute transmissivity t
    nlayers = size(t, 2)
    NF = eltype(t)

    (; dry_absorption, water_vapor_absorption, co2_sensitivity, reference_co2) = transmissivity
    (; reference_pressure) = model.atmosphere
    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]
    has_humidity = haskey(vars.grid, :humidity)     # dry models: only the dry optical depth

    # dry optical depth per p₀, logarithmic in CO₂, clamped to be non-negative for very low CO₂
    co2 = co2_concentration(vars.prognostic, reference_co2)
    dry = max(dry_absorption + co2_sensitivity * log(co2 / reference_co2), zero(NF))

    for k in 1:nlayers
        q = has_humidity ?
            max(get_prognostic_step(vars.grid.humidity, model.time_stepping, transmissivity)[ij, k], zero(NF)) :
            zero(NF)
        Δp = pressure_thickness(k, pₛ, coord)
        dτ = (dry + water_vapor_absorption * q) * Δp / reference_pressure
        t[ij, k] = exp(-dτ)
    end
    return t
end
