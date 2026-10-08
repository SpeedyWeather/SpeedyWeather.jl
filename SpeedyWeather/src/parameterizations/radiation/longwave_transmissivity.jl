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

export CloudyLongwaveTransmissivity

"""Longwave transmissivity of a clear-sky transmissivity `clear_sky` and the clouds of a
prognostic cloud scheme, e.g. [`PrognosticCloudCondensation`](@ref). In every layer the clear-sky
transmissivity is multiplied by `1 - C ε` with the cloud fraction `C` and the emissivity of the
cloudy part

    ε = 1 - exp(-D (κₗ Wₗ + κᵢ Wᵢ))

from the in-cloud liquid and ice water paths `W` [kg/m²], the diffusivity factor `D` and the mass
absorption coefficients `κₗ` (constant) and `κᵢ = a + b/rᵢ` with the ice effective radius `rᵢ`
(Kiehl et al. 1998, CCM3, after Ebert and Curry 1992). Clouds of different layers are
independent (random overlap). Without a prognostic cloud scheme the cloud state is zero and
this is the clear-sky transmissivity. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct CloudyLongwaveTransmissivity{NF, T} <: AbstractLongwaveTransmissivity
    "[OPTION] Clear-sky longwave transmissivity"
    @component clear_sky::T

    "[OPTION] Diffusivity factor of the longwave emissivity [1]"
    @param diffusivity::NF = 1.66 (bounds = Positive,)

    "[OPTION] Mass absorption coefficient of cloud liquid [m²/kg]"
    @param liquid_mass_absorption::NF = 90.361 (bounds = Nonnegative,)

    "[OPTION] Mass absorption coefficient of cloud ice, constant part a [m²/kg]"
    @param ice_mass_absorption::NF = 5 (bounds = Nonnegative,)

    "[OPTION] Mass absorption coefficient of cloud ice, part b/rᵢ inverse to the effective radius [m³/kg]"
    @param ice_mass_absorption_radius::NF = 1.0e-3 (bounds = Nonnegative,)
end

Adapt.@adapt_structure CloudyLongwaveTransmissivity

function CloudyLongwaveTransmissivity(SG::SpectralGrid; clear_sky = FriersonLongwaveTransmissivity(SG), kwargs...)
    return CloudyLongwaveTransmissivity{SG.NF, typeof(clear_sky)}(; clear_sky, kwargs...)
end

initialize!(transmissivity::CloudyLongwaveTransmissivity, model::AbstractModel) =
    initialize!(transmissivity.clear_sky, model)

variables(transmissivity::CloudyLongwaveTransmissivity) = (
    cloud_state_variables()...,
    variables(transmissivity.clear_sky)...,
    ParameterizationVariable(:outgoing_longwave_clear_sky, Grid2D(), desc = "TOA longwave radiation up without clouds", units = "W/m^2"),
)

@propagate_inbounds function transmissivity!(ij, vars, transmissivity::CloudyLongwaveTransmissivity, model)
    # clear-sky transmissivity first, into the scratch array it returns, a copy in scratch b for
    # the clear-sky outgoing longwave (see clear_sky_outgoing_longwave!)
    t = transmissivity!(ij, vars, transmissivity.clear_sky, model)
    t_clear = vars.scratch.grid.b

    (; cloud_fraction, cloud_liquid_water, cloud_ice_water, cloud_ice_effective_radius) = vars.parameterizations
    (; diffusivity, liquid_mass_absorption, ice_mass_absorption, ice_mass_absorption_radius) = transmissivity
    NF = eltype(t)
    nlayers = size(t, 2)

    pₛ = vars.parameterizations.surface_pressure[ij]
    coord = model.geometry.vertical_coordinates
    g = model.planet.gravity

    for k in 1:nlayers
        C = cloud_fraction[ij, k]
        in_cloud_mass = pressure_thickness(k, pₛ, coord) / g / max(C, eps(NF))  # [kg/m²] per kg/kg in the cloudy part
        κᵢ = ice_mass_absorption + ice_mass_absorption_radius / max(cloud_ice_effective_radius[ij, k], NF(1.0e-6))
        absorption = liquid_mass_absorption * cloud_liquid_water[ij, k] + κᵢ * cloud_ice_water[ij, k]
        emissivity = 1 - exp(-diffusivity * absorption * in_cloud_mass)
        t_clear[ij, k] = t[ij, k]
        t[ij, k] *= 1 - C * emissivity
    end
    return t
end

# with cloudy transmissivity also the clear-sky outgoing longwave from the clear-sky transmissivity
# that CloudyLongwaveTransmissivity keeps in scratch b
@propagate_inbounds function parameterization!(ij, vars, radiation::OneBandLongwave{<:CloudyLongwaveTransmissivity}, model)
    t = transmissivity!(ij, vars, radiation.transmissivity, model)
    longwave_radiative_transfer!(ij, vars, t, radiation.radiative_transfer, model)
    clear_sky_outgoing_longwave!(ij, vars, vars.scratch.grid.b, radiation.radiative_transfer, model)
    return nothing
end
