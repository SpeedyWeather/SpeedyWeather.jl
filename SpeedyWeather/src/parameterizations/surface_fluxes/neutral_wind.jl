abstract type AbstractNeutralWindSpeed <: AbstractParameterization end

export NeutralWindSpeed

"""Ocean-based neutral wind speed calculation from actual wind speed, 
derived from ERA5 data via symbolic regression. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NeutralWindSpeed{NF} <: AbstractNeutralWindSpeed
    "[OPTION] symbolic regression parameter 1 [1]"
    @param c1::NF = -0.039317116

    "[OPTION] symbolic regression parameter 2 [1]"
    @param c2::NF = -2.9858496

    "[OPTION] symbolic regression parameter 3 [1]"
    @param c3::NF = 2.0046231e-10

    "[OPTION] symbolic regression parameter 4 [1]"
    @param c4::NF = 1.0768474

    "[OPTION] symbolic regression parameter 5 [1]"
    @param c5::NF = 0.20268184

    "[OPTION] symbolic regression parameter 6 [1]"
    @param c6::NF = 1.2684147

    "[OPTION] symbolic regression parameter 7 [1]"
    @param c7::NF = -0.94933933

    "[OPTION] symbolic regression parameter 8 [1]"
    @param c8::NF = 0.041551278

    "[OPTION] symbolic regression parameter 9 [1]"
    @param c9::NF = 5.8649142

    # Safety parameters
    "[OPTION] Minimum wind speed to avoid zero surface fluxes [m/s]"
    @param minimum_wind::NF = 1.0e-6

    "[OPTION] Minimum argument for logarithm to avoid log(0) [1]"
    @param log_safety::NF = 1.0e-8
end

# every neutral wind speed parameterization is expected to write this on every time step
variables(::AbstractNeutralWindSpeed) = (
    ParameterizationVariable(:neutral_wind_speed, GridXY(), desc = "Neutral surface wind speed", units = "m/s"),
)

Adapt.@adapt_structure NeutralWindSpeed
NeutralWindSpeed(SG::SpectralGrid; kwargs...) = NeutralWindSpeed{SG.NF}(; kwargs...)
initialize!(::NeutralWindSpeed, ::PrimitiveEquation) = nothing

@propagate_inbounds function parameterization!(ij, vars, nw::NeutralWindSpeed{NF}, model) where {NF}
    (; land_fraction) = model.land_sea_mask
    (land_fraction[ij] < 1) || return nothing # TODO train a land-based neutral wind speed parameterization
    return neutral_wind_speed(ij, vars, nw, model)
end

@propagate_inbounds function neutral_wind_speed(ij, vars, nw::NeutralWindSpeed{NF}, model) where {NF}
    (; surface_wind_speed) = vars.parameterizations
    (; surface_air_temperature) = vars.parameterizations

    sst = vars.prognostic.ocean.sea_surface_temperature[ij]
    t_diff = surface_air_temperature[ij] - sst # TODO: replace SST with ocean skin temperature
    ws_safe = max(surface_wind_speed[ij], nw.minimum_wind)
    log_arg = max(nw.c1 * ws_safe * t_diff + exp(t_diff), nw.log_safety)

    numerator = 2 * t_diff + nw.c8 * exp(t_diff) - nw.c3 * (nw.c4^surface_air_temperature[ij])
    denominator = t_diff * (log(log_arg) + nw.c2) + nw.c5 * (nw.c6^ws_safe) + nw.c9 * (ws_safe^nw.c7) + ws_safe

    vars.parameterizations.neutral_wind_speed[ij] = max(surface_wind_speed[ij] - (numerator / denominator), 0)

    return nothing
end

export NonNeutralSurfaceWindSpeed

"""Approximates the neutral wind speed by the (non-neutral) surface wind speed,
i.e. no stability correction. Fallback for parameterizations that require
a neutral wind speed, e.g. `LearnedOceanRoughness`. Does not compute or store anything,
`get_neutral_wind_speed` returns the surface wind speed instead. Fields are $(TYPEDFIELDS)"""
struct NonNeutralSurfaceWindSpeed <: AbstractNeutralWindSpeed end

NonNeutralSurfaceWindSpeed(::SpectralGrid) = NonNeutralSurfaceWindSpeed()
initialize!(::NonNeutralSurfaceWindSpeed, ::PrimitiveEquation) = nothing
variables(::NonNeutralSurfaceWindSpeed) = ()        # no neutral wind speed array needed
parameterization!(ij, vars, ::NonNeutralSurfaceWindSpeed, model) = nothing

"""$(TYPEDSIGNATURES)
Neutral wind speed at grid point `ij` as computed by the neutral wind speed scheme,
or approximated by the surface wind speed for `NonNeutralSurfaceWindSpeed`."""
@propagate_inbounds get_neutral_wind_speed(ij, vars, ::AbstractNeutralWindSpeed) =
    vars.parameterizations.neutral_wind_speed[ij]
@propagate_inbounds get_neutral_wind_speed(ij, vars, ::NonNeutralSurfaceWindSpeed) =
    vars.parameterizations.surface_wind_speed[ij]
