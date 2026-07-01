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
    pₛ = vars.parameterizations.surface_pressure[ij]    # actual surface pressure, see NOTE below

    for k in 1:nlayers
        Δσₖ = pressure_thickness(k, pₛ, coord) / pₛ   # sigma thickness of layer k
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
    pₛ = vars.parameterizations.surface_pressure[ij]    # actual surface pressure, see NOTE below
    j = model.geometry.whichring[ij]
    sinlat = model.geometry.sinlat[j]   # sin(latitude), avoids calling sind() on traced scalars

    # Frierson 2006, eq. (4), (5) but in a differential form, computing dτ between half levels below and above
    # --- τ(k=1/2)                  # half level above
    # dt = τ(k=1+1/2) - τ(k=1/2)    # differential optical depth on layer k
    # --- τ(k=1+1/2)                # half level below

    τ_above::NF = 0
    τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2
    for k in 1:nlayers              # loop over half levels below
        σₖ = pressure_below(k, pₛ, coord) / pₛ  # sigma at half level below full level k
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
transmissivity bit-identical across CPU/GPU/Reactant. Only valid for state-independent
coordinates (`SigmaCoordinates`); state-dependent coordinates compute `transmissivity!` per
column at run time instead."""
function fill_longwave_transmissivity!(t::AbstractMatrix, CLT::ConstantLongwaveTransmissivity, model)
    nlayers, nlat = size(t)
    # pull the coordinate to host so the loop runs host libm `exp`, not device/XLA `exp`;
    # on CPU this is a no-op, under Reactant it transfers ConcretePJRTArray → Vector
    coord = on_architecture(CPU(), model.geometry.vertical_coordinates)
    pₛ = model.atmosphere.reference_pressure # reference (not prognostic) surface pressure, see NOTE above
    τ = -log(CLT.transmissivity)            # total optical depth of the atmosphere
    for j in 1:nlat, k in 1:nlayers
        Δσₖ = pressure_thickness(k, pₛ, coord) / pₛ   # sigma thickness of layer k at the reference pₛ
        t[k, j] = exp(-τ * Δσₖ)             # transmissivity through layer k (lat-independent)
    end
    return t
end

function fill_longwave_transmissivity!(t::AbstractMatrix, transmissivity::FriersonLongwaveTransmissivity, model)
    nlayers, nlat = size(t)
    NF = eltype(t)
    (; τ₀_equator, τ₀_pole, fₗ) = transmissivity
    # pull the coordinate + sinlat to host so the loop runs host libm `exp`, not device/XLA `exp`
    coord = on_architecture(CPU(), model.geometry.vertical_coordinates)
    pₛ = model.atmosphere.reference_pressure # reference (not prognostic) surface pressure, see NOTE above
    sinlat = on_architecture(CPU(), model.geometry.sinlat)          # sin(latitude) per ring
    for j in 1:nlat
        τ_above = zero(NF)
        τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat[j]^2
        for k in 1:nlayers                  # loop over half levels below
            σₖ = pressure_below(k, pₛ, coord) / pₛ  # sigma at half level below full level k at the reference pₛ
            τ_below = τ₀ * (fₗ * σₖ + (1 - fₗ) * σₖ^4)
            t[k, j] = exp(-(τ_below - τ_above))
            τ_above = τ_below
        end
    end
    return t
end

# Provide the per-column layer transmissivity for grid point `ij` in the scratch array (indexed
# `[ij, k]`), which the radiative transfer then reads. Dispatches on the vertical coordinate:
#  - `SigmaCoordinates` (state-independent): copy the precomputed per-latitude-ring column
#  - any other coordinate (state-dependent): compute on the fly from the actual surface pressure
@propagate_inbounds function transmissivity_column!(ij, vars, transmissivity, ::SigmaCoordinates, model)
    t = vars.scratch.grid.a
    t_precomputed = vars.parameterizations.longwave_transmissivity   # nlayers × nlat
    j = model.geometry.whichring[ij]
    nlayers = size(t, 2)
    for k in 1:nlayers
        t[ij, k] = t_precomputed[k, j]
    end
    return t
end

@propagate_inbounds transmissivity_column!(ij, vars, transmissivity, ::AbstractVerticalCoordinates, model) =
    transmissivity!(ij, vars, transmissivity, model)
