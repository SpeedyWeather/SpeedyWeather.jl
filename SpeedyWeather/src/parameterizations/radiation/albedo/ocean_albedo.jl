# OceanSeaIceAlbedo
export OceanSeaIceAlbedo

"""Albedo that scales linearly between ocean and ice albedo depending on sea ice concentration.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OceanSeaIceAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Albedo over open ocean [1]"
    @param albedo_ocean::NF = 0.06 (bounds = 0 .. 1,)

    "[OPTION] Albedo over sea ice at concentration=1 [1]"
    @param albedo_ice::NF = 0.6 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure OceanSeaIceAlbedo
OceanSeaIceAlbedo(SG::SpectralGrid; kwargs...) = OceanSeaIceAlbedo{SG.NF}(; kwargs...)
initialize!(::OceanSeaIceAlbedo, ::PrimitiveEquation) = nothing

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::OceanSeaIceAlbedo, model)
    (; albedo_ocean, albedo_ice) = scheme
    NF = eltype(albedo)
    ℵ = haskey(vars.prognostic.ocean, :sea_ice_concentration) ? vars.prognostic.ocean.sea_ice_concentration[ij] : zero(NF)

    # set ocean albedo linearly between ocean and ice depending on sea ice concentration
    return albedo[ij] = albedo_ocean + ℵ * (albedo_ice - albedo_ocean)
end

"""Whitecap (foam) correction of the ocean albedo, a component of `JinOceanAlbedo` (not an albedo
itself). Subtypes `F <: AbstractFoam` are callable as `F(α, w)` and return the albedo `α` [1]
corrected for whitecaps at wind speed `w` [m/s]."""
abstract type AbstractFoam <: AbstractModelComponent end

# components of JinOceanAlbedo need no initialization unless a method is defined for them
initialize!(::AbstractFoam, ::PrimitiveEquation) = nothing

export WhitecapFoam, NoFoam

"""Whitecap (foam) correction following Koepke (1984) as used by Jin et al. (2011) eqs. (16)–(17):
constant foam albedo and fractional whitecap coverage f = a * w^b from wind speed w [m/s].
Component of `JinOceanAlbedo`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct WhitecapFoam{NF} <: AbstractFoam
    "[OPTION] Albedo of whitecaps/foam [1], Koepke (1984)"
    @param albedo_foam::NF = 0.55 (bounds = 0 .. 1,)

    "[OPTION] Whitecap coverage coefficient [(s/m)^b], Jin et al. eq. (16)"
    @param coverage_coefficient::NF = 2.95e-6 (bounds = Positive,)

    "[OPTION] Whitecap coverage wind speed exponent [1], Jin et al. eq. (16)"
    @param coverage_exponent::NF = 3.52 (bounds = Positive,)
end

Adapt.@adapt_structure WhitecapFoam
WhitecapFoam(SG::SpectralGrid; kwargs...) = WhitecapFoam{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES) Fractional whitecap coverage in [0, 1] at wind speed `w` [m/s], Jin et al. (2011), eq. (16)."""
@inline function whitecap_coverage(foam::WhitecapFoam, w)
    (; coverage_coefficient, coverage_exponent) = foam
    return min(coverage_coefficient * w^coverage_exponent, one(w))
end

"""$(TYPEDSIGNATURES) Albedo `α` corrected for whitecaps at wind speed `w` [m/s], Jin et al. (2011), eq. (17)."""
@inline function (foam::WhitecapFoam)(α, w)
    f = whitecap_coverage(foam, w)
    return f * foam.albedo_foam + (1 - f) * α
end

# Correct both albedos together so WhitecapFoam computes coverage only once.
@inline foam_albedos(foam::AbstractFoam, black_sky, white_sky, w) = (foam(black_sky, w), foam(white_sky, w))

@inline function foam_albedos(foam::WhitecapFoam, black_sky, white_sky, w)
    f = whitecap_coverage(foam, w)
    black_sky = f * foam.albedo_foam + (1 - f) * black_sky
    white_sky = f * foam.albedo_foam + (1 - f) * white_sky
    return black_sky, white_sky
end

"""No whitecap correction, component of `JinOceanAlbedo`."""
struct NoFoam <: AbstractFoam end
Adapt.@adapt_structure NoFoam
NoFoam(SG::SpectralGrid) = NoFoam()

"""$(TYPEDSIGNATURES) Albedo `α` unchanged, no whitecap correction depending on wind speed `w`."""
@inline (::NoFoam)(α, w) = α

# JinOceanAlbedo is composed of components with three different roles: surface scattering,
# ocean volume scattering and whitecaps (foam). Each role has its own abstract type (none of
# them are albedos themselves, i.e. they cannot be used as `ocean` in `OceanLandAlbedo`)
# so that the components of `JinOceanAlbedo` can be constrained to their role.

export JinOceanAlbedo

"""Ocean albedo following Jin et al. (2011, Optics Express 19(27), 26429), depending on the
solar zenith angle and the surface wind speed. The direct (black-sky) and diffuse (white-sky)
albedos are computed from the surface (Fresnel) reflection of the air-water interface,
roughened by wind, plus the volume scattering of the ocean below, with a correction for
whitecaps (foam) at high wind speeds. The shortwave radiation scheme blends both with the
fraction of direct radiation. Sea ice is blended in linearly with its concentration
as in `OceanSeaIceAlbedo`. To be used as `ocean` albedo in `OceanLandAlbedo` (the default).
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct JinOceanAlbedo{
    NF,
    SDIR <: AbstractOceanSurfaceScattering,
    SDIF <: AbstractOceanSurfaceScattering,
    V <: AbstractOceanVolumeScattering,
    F <: AbstractFoam
} <: AbstractAlbedo

    "Ratio of air-water refractive indices"
    @param n::NF = 1.34

    "Refractive index of water in visible spectrum"
    @param n₀::NF = 1.34

    "Effective cosine for diffuse incidence, Jin et al. eq. (10)"
    @param μₑ::NF = 0.676

    "[OPTION] Albedo over sea ice at concentration=1 [1], same default as in `OceanSeaIceAlbedo`"
    @param albedo_ice::NF = 0.6 (bounds = 0 .. 1,)

    "Direct-beam surface reflection"
    @component ocean_direct_reflection::SDIR

    "Diffuse surface reflection"
    @component ocean_diffuse_reflection::SDIF

    "[OPTION] Ocean volume scattering, `ConstantOceanVolumeScattering` (`SpectralOceanVolumeScattering` is future work)"
    @component ocean_volume::V

    "[OPTION] Whitecap (foam) correction, `WhitecapFoam` or `NoFoam`"
    @component foam::F
end

Adapt.@adapt_structure JinOceanAlbedo
function JinOceanAlbedo(
        SG::SpectralGrid;
        ocean_direct_reflection = OceanDirectReflection(SG),
        ocean_diffuse_reflection = OceanDiffuseReflection(SG),
        ocean_volume = ConstantOceanVolumeScattering(SG),
        foam = WhitecapFoam(SG),
        kwargs...
    )
    # construct first so that the type constraints on ocean_volume and foam are checked first
    albedo = JinOceanAlbedo{SG.NF, typeof(ocean_direct_reflection), typeof(ocean_diffuse_reflection), typeof(ocean_volume), typeof(foam)}(;
        ocean_direct_reflection, ocean_diffuse_reflection, ocean_volume, foam, kwargs...
    )
    check_ocean_volume_supported(albedo.ocean_volume)
    return albedo
end

# SpectralOceanVolumeScattering is future work: no shortwave scheme computes a surface albedo
# per band yet. Checked in the constructor and in initialize! to also catch the @kwdef constructor
check_ocean_volume_supported(::AbstractOceanVolumeScattering) = nothing
check_ocean_volume_supported(::SpectralOceanVolumeScattering) = throw(
    ArgumentError(
        "SpectralOceanVolumeScattering is future work and cannot be used yet: it requires a spectral " *
            "shortwave radiation scheme with a surface albedo per band. Use ConstantOceanVolumeScattering instead."
    )
)

function variables(::JinOceanAlbedo)
    return (
        ParameterizationVariable(:black_sky_albedo, GridXY(), desc = "Black-sky ocean albedo", units = "1", namespace = :ocean),
        ParameterizationVariable(:white_sky_albedo, GridXY(), desc = "White-sky ocean albedo", units = "1", namespace = :ocean),
    )
end

function initialize!(albedo::JinOceanAlbedo, model::PrimitiveEquation)
    check_ocean_volume_supported(albedo.ocean_volume)
    initialize!(albedo.ocean_direct_reflection, model)
    initialize!(albedo.ocean_diffuse_reflection, model)
    initialize!(albedo.ocean_volume, model)
    initialize!(albedo.foam, model)
    return nothing
end

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::JinOceanAlbedo, model)
    (; ocean) = vars.parameterizations
    (; n, n₀, μₑ) = scheme
    NF = eltype(albedo)

    # Cosine of solar zenith angle, only when the sun is up
    μ = clamp(vars.parameterizations.cos_zenith[ij], zero(NF), one(NF))
    w = vars.parameterizations.surface_wind_speed[ij]
    σ = sqrt(NF(0.003) + NF(0.00512) * w)     # Cox & Munk, Jin et al. eq. (2)

    p = (; μ, σ, n, n₀)

    # Surface (Fresnel) reflection, Jin eqs. (1) and (5b)
    αs_dir = scheme.ocean_direct_reflection(p)
    αs_dif = scheme.ocean_diffuse_reflection(p)

    # Ocean volume scattering, Jin eqs. (6) and (10)
    αw_dir = scheme.ocean_volume((; p..., αₛ = αs_dir), μ)
    αw_dif = scheme.ocean_volume((; p..., αₛ = αs_dif), μₑ)

    # Direct (black-sky) and diffuse (white-sky) albedo, Jin eqs. (11)–(12),
    # blended by the shortwave scheme via combine_surface_albedo!
    black_sky = αs_dir + αw_dir
    white_sky = αs_dif + αw_dif

    # Whitecap (foam) correction, Jin et al. eqs. (16)–(17)
    black_sky, white_sky = foam_albedos(scheme.foam, black_sky, white_sky, w)

    ℵ = haskey(vars.prognostic.ocean, :sea_ice_concentration) ? vars.prognostic.ocean.sea_ice_concentration[ij] : zero(NF)
    black_sky += ℵ * (scheme.albedo_ice - black_sky)
    white_sky += ℵ * (scheme.albedo_ice - white_sky)

    # Without a direct beam (sun below the horizon) the black-sky albedo would be the Fresnel
    # reflection at grazing incidence (0.3-0.6), use the diffuse white-sky albedo instead.
    # No effect on the shortwave fluxes (zero at night), but keeps the albedo output meaningful
    black_sky = ifelse(μ > 0, black_sky, white_sky)

    ocean.black_sky_albedo[ij] = black_sky
    ocean.white_sky_albedo[ij] = white_sky

    # Fallback for radiation schemes with no direct/diffuse radiation fraction
    albedo[ij] = white_sky

    return nothing
end