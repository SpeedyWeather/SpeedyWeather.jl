abstract type AbstractOceanScattering <: AbstractScattering end

"""Surface reflection of the ocean, a component of `JinOceanAlbedo` (not an albedo itself).
Subtypes `R <: AbstractOceanReflection` are callable as `R(p)` with `p` a NamedTuple of the
cosine of the solar zenith angle `μ`, the surface roughness `σ` and the refractive indices
`n`, `n₀` and return the surface reflection albedo [1]."""
abstract type AbstractOceanSurfaceScattering <: AbstractOceanScattering end

"""Volume scattering of the ocean (light scattered back from below the surface), a component of
`JinOceanAlbedo` (not an albedo itself). Subtypes `V <: AbstractOceanVolumeScattering` are callable as
`V(p, μ)` with `p` a NamedTuple of the surface roughness `σ` and the surface reflection albedo `αₛ`,
and `μ` the cosine of the (effective) incidence angle, and return the volume albedo [1]."""
abstract type AbstractOceanVolumeScattering <: AbstractOceanScattering end

# components of JinOceanAlbedo need no initialization unless a method is defined for them
initialize!(::AbstractOceanSurfaceScattering, ::PrimitiveEquation) = nothing
initialize!(::AbstractOceanVolumeScattering, ::PrimitiveEquation) = nothing

"""Direct-beam surface reflection of the ocean following Jin et al. (2011), eq. (1): the
Fresnel reflectance of the air-water interface with a correction for the roughness of the
wind-driven surface, a regression in the cosine of the solar zenith angle μ and the surface
roughness σ. Component of `JinOceanAlbedo`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OceanDirectReflection{NF} <: AbstractOceanSurfaceScattering
    "Regression coefficient from Jin et al."
    @param c0::NF = 0.0152

    "Regression coefficient from Jin et al."
    @param c1::NF = -1.7873

    "Regression coefficient from Jin et al."
    @param c2::NF = 6.8972

    "Regression coefficient from Jin et al."
    @param c3::NF = -8.5778

    "Regression coefficient from Jin et al."
    @param c4::NF = 4.071

    "Regression coefficient from Jin et al."
    @param c5::NF = -7.6446

    "Regression coefficient from Jin et al."
    @param c6::NF = 0.1643

    "Regression coefficient from Jin et al."
    @param c7::NF = -7.8409

    "Regression coefficient from Jin et al."
    @param c8::NF = -3.5639

    "Regression coefficient from Jin et al."
    @param c9::NF = -2.3588

    "Regression coefficient from Jin et al."
    @param c10::NF = 10.0538
end

Adapt.@adapt_structure OceanDirectReflection
OceanDirectReflection(SG::SpectralGrid; kwargs...) = OceanDirectReflection{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Unpolarized Fresnel reflectance of a flat air–water interface for light incident from air
with cosine of incidence angle `μ` onto a medium of relative refractive index `n > 1`.
Average of the perpendicular (s) and parallel (p) polarized reflectances. Returns
((n-1)/(n+1))² at normal incidence (μ = 1) and 1 at grazing incidence (μ = 0)."""
@inline function fresnel_reflectance(n, μ)
    s = sqrt(n^2 - 1 + μ^2)                 # n cos θₜ via Snell's law, θₜ refraction angle
    r_s = (μ - s) / (μ + s)                 # perpendicular (s) polarization
    r_p = (n^2 * μ - s) / (n^2 * μ + s)     # parallel (p) polarization
    return (r_s^2 + r_p^2) / 2
end

"""$(TYPEDSIGNATURES)
Direct-beam surface reflection albedo [1] following Jin et al. (2011), eq. (1),
`r(n, μ) - r(n, μ)/r(n₀, μ) * f(μ, σ)` with `r` the Fresnel reflectance and `f(μ, σ)`
the regression that corrects for the roughness of a windy, wavy surface. The cosine of the
solar zenith angle `μ`, the surface roughness `σ` and the refractive indices `n`, `n₀`
are taken from the NamedTuple `p`."""
@inline function (c::OceanDirectReflection)(p)
    (; μ, σ, n, n₀) = p
    (; c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) = c

    # Roughness correction f(μ, σ)
    term1 = c0 + c1 * μ + c2 * μ^2 + c3 * μ^3 + c4 * σ + c5 * μ * σ
    term2 = exp(c6 + c7 * μ + c8 * μ^2 + c9 * σ + c10 * μ * σ)
    f = term1 * term2

    # Now fresnel reflectances
    r = fresnel_reflectance(n, μ)
    r₀ = fresnel_reflectance(n₀, μ)

    albedo_direct_reflection = r - (r / r₀) * f

    return albedo_direct_reflection
end

"""Diffuse surface reflection of the ocean following Jin et al. (2011), eq. (5b), a linear
regression in the refractive index n and the surface roughness σ. This is the parameterization
for diffuse radiation under clouds rather than for isotropic diffuse radiation, as SpeedyWeather
splits the shortwave radiation into direct and diffuse by cloud cover.
Component of `JinOceanAlbedo`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OceanDiffuseReflection{NF} <: AbstractOceanSurfaceScattering
    "Regression coefficient from Jin et al."
    @param c0::NF = -0.1479

    "Regression coefficient from Jin et al."
    @param c1::NF = 0.1502

    "Regression coefficient from Jin et al."
    @param c2::NF = -0.0176
end

Adapt.@adapt_structure OceanDiffuseReflection
OceanDiffuseReflection(SG::SpectralGrid; kwargs...) = OceanDiffuseReflection{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Diffuse surface reflection albedo [1] following Jin et al. (2011), eq. (5b), from the
refractive index `n` and the surface roughness `σ` in the NamedTuple `p`."""
@inline function (c::OceanDiffuseReflection)(p)
    (; σ, n) = p
    (; c0, c1, c2) = c
    return c0 + c1 * n + c2 * n * σ
end

export ConstantOceanVolumeScattering

"""Constant broadband volume-scattering albedo of the ocean for case 1 waters, Jin et al. (2011),
eq. (15), independent of the incidence angle. For shortwave schemes without spectral bands.
Component of `JinOceanAlbedo`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantOceanVolumeScattering{NF} <: AbstractOceanVolumeScattering
    "[OPTION] Broadband volume-scattering albedo for case 1 waters, Jin et al. eq. (15) [1]"
    @param albedo::NF = 0.006 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure ConstantOceanVolumeScattering
ConstantOceanVolumeScattering(SG::SpectralGrid; kwargs...) = ConstantOceanVolumeScattering{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES) Constant volume-scattering albedo [1], independent of `p` and `μ`."""
@inline (c::ConstantOceanVolumeScattering)(p, μ) = c.albedo

"""Spectral volume-scattering albedo of the ocean following Jin et al. (2011), eqs. (6)–(9),
with the backscattering and absorption of case 1 waters from Morel & Maritorena (2001)
depending on the chlorophyll concentration, evaluated at the centre `wavelengths` [nm]
of the shortwave bands.

!!! warning "Future work, not yet usable"
    This requires a spectral shortwave radiation scheme that computes a surface albedo per band
    and passes the band index as `band` in the NamedTuple `p` when calling `V(p, μ)`. No shortwave
    scheme in SpeedyWeather does this yet, so `JinOceanAlbedo` throws an `ArgumentError` when
    constructed or initialized with a `SpectralOceanVolumeScattering`.
    Use `ConstantOceanVolumeScattering` instead.

Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct SpectralOceanVolumeScattering{NF, VectorType} <: AbstractOceanVolumeScattering
    "[OPTION] Chlorophyll concentration [mg/m³], clamped to ≥ 0.02 (M&M eq. 14 lower limit)"
    @param chl::NF = 0.1 (bounds = Positive,)

    "Regression coefficient from Jin et al."
    @param r0::NF = 0.4817

    "Regression coefficient from Jin et al."
    @param r1::NF = 0.0149

    "Regression coefficient from Jin et al."
    @param r2::NF = 0.207

    "Regression coefficient from Jin et al."
    @param β0::NF = 0.6279

    "Regression coefficient from Jin et al."
    @param β1::NF = 0.2227

    "Regression coefficient from Jin et al."
    @param β2::NF = 0.0513

    "Regression coefficient from Jin et al."
    @param β3::NF = 0.2465

    "Regression coefficient from Jin et al."
    @param β4::NF = 0.3119

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b0::NF = 0.002

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b1::NF = 0.01

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b2::NF = 0.5

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b3::NF = 0.25

    "Coefficient from Morel & Maritorena, eqn (13), reference wavelength [nm]"
    @param b4::NF = 550

    "Coefficient from Morel & Maritorena, eqn (12)"
    @param b5::NF = 0.416

    "Coefficient from Morel & Maritorena, eqn (12)"
    @param b6::NF = 0.766

    "Coefficient from Morel & Maritorena, eqn (14)"
    @param b7::NF = 0.3

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a0::NF = 0.06

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a1::NF = 0.65

    "Coefficient from Morel & Maritorena, eqn (17)"
    @param a2::NF = 0.014

    "Coefficient from Morel & Maritorena, eqn (17), reference wavelength [nm]"
    @param a3::NF = 440

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a4::NF = 0.2

    "Pure water absorption at 440 nm in the CDOM term, Pope & Fry (1997) as in Jin et al.'s code [m⁻¹]"
    @param a_w_440::NF = 0.00635

    "Band centre wavelengths [nm]"
    wavelengths::VectorType

    "Per-band R₀ intercept A in R₀ = A + Bμ, computed in initialize!"
    R₀_intercept::VectorType

    "Per-band R₀ slope B in R₀ = A + Bμ, computed in initialize!"
    R₀_slope::VectorType
end

Adapt.@adapt_structure SpectralOceanVolumeScattering
function SpectralOceanVolumeScattering(SG::SpectralGrid; wavelengths = [450, 550, 650], kwargs...)
    (; NF, VectorType, architecture) = SG
    nbands = length(wavelengths)
    return SpectralOceanVolumeScattering{NF, VectorType}(;
        wavelengths = on_architecture(architecture, NF.(wavelengths)),
        R₀_intercept = on_architecture(architecture, zeros(NF, nbands)),
        R₀_slope = on_architecture(architecture, zeros(NF, nbands)),
        kwargs...,
    )
end

"""$(TYPEDSIGNATURES) Linear interpolation of tabulated `ys` at ascending `xs`,
clamped to the end values outside the table (no extrapolation)."""
@inline function interpolate_table(xs::NTuple{N}, ys::NTuple{N}, x) where {N}
    x <= xs[1] && return oftype(x, ys[1])
    x >= xs[N] && return oftype(x, ys[N])
    lo, hi = 1, N                       # binary search for xs[lo] ≤ x < xs[hi]
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        xs[mid] <= x ? (lo = mid) : (hi = mid)
    end
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return oftype(x, (1 - t) * ys[lo] + t * ys[hi])
end

# Normalized phytoplankton absorption A_chl(λ), A_chl(440) = 1, from Jin et al. (2011) supplementary code
const ACHL_WAVELENGTHS = Tuple(350.0:10.0:700.0)
const ACHL = (
    0.56, 0.57, 0.59, 0.61, 0.66, 0.687, 0.828, 0.913, 0.973, 1.0, 0.944, 0.917,
    0.87, 0.798, 0.75, 0.668, 0.618, 0.528, 0.474, 0.416, 0.357, 0.294, 0.276, 0.291,
    0.282, 0.236, 0.252, 0.276, 0.317, 0.334, 0.356, 0.441, 0.595, 0.502, 0.329, 0.215,
)

# Pure water absorption a_w(λ) [m⁻¹], from Jin et al. (2011) supplementary code:
# Smith & Baker (1981) for λ ≤ 510 nm, Pope & Fry (1997) for 512.5–727.5 nm
const AW_WAVELENGTHS = (Tuple(200.0:10.0:510.0)..., Tuple(512.5:2.5:730.0)..., Tuple(740.0:10.0:800.0)...)
const AW = (
    3.07, 1.99, 1.31, 0.927, 0.72, 0.559, 0.457, 0.373,
    0.288, 0.215, 0.141, 0.105, 0.0844, 0.0678, 0.0561, 0.0463,
    0.0379, 0.03, 0.022, 0.0191, 0.0171, 0.0162, 0.0153, 0.0144,
    0.0145, 0.0145, 0.0156, 0.0156, 0.0176, 0.0196, 0.0257, 0.0357,
    0.0372, 0.0396, 0.0399, 0.0409, 0.0416, 0.0417, 0.0428, 0.0434,
    0.0447, 0.0452, 0.0466, 0.0474, 0.0489, 0.0511, 0.0537, 0.0565,
    0.0593, 0.0596, 0.0606, 0.0619, 0.064, 0.0642, 0.0672, 0.0695,
    0.0733, 0.0772, 0.0836, 0.0896, 0.0989, 0.11, 0.122, 0.1351,
    0.1516, 0.1672, 0.1925, 0.2224, 0.247, 0.2577, 0.2629, 0.2644,
    0.2665, 0.2678, 0.2707, 0.2755, 0.281, 0.2834, 0.2904, 0.2916,
    0.2995, 0.3012, 0.3077, 0.3108, 0.322, 0.325, 0.335, 0.34,
    0.358, 0.371, 0.393, 0.41, 0.424, 0.429, 0.436, 0.439,
    0.448, 0.448, 0.461, 0.465, 0.478, 0.486, 0.502, 0.516,
    0.538, 0.559, 0.592, 0.624, 0.663, 0.704, 0.756, 0.827,
    0.914, 1.007, 1.119, 1.231, 1.356, 1.489, 1.678, 1.799,
    2.38, 2.47, 2.55, 2.51, 2.36, 2.16, 2.07,
)

@inline A_chl(λnm) = interpolate_table(ACHL_WAVELENGTHS, ACHL, λnm)
@inline a_w(λnm) = interpolate_table(AW_WAVELENGTHS, AW, λnm)

# From Morel (1974) and Bricaud (1998)
@inline b_w(λnm) = oftype(λnm, 0.00288) * (λnm / 500)^oftype(λnm, -4.32)

"""$(TYPEDSIGNATURES) Coefficients A, B of the subsurface irradiance reflectance R₀ = A + Bμ
at wavelength `λ` [nm] for chlorophyll `chl` [mg/m³], Jin et al. eqs. (8)–(9) with
Morel & Maritorena (2001) eqs. (11)–(18). Called once per band at initialization."""
function R₀_coefficients(c::SpectralOceanVolumeScattering, λ, chl)
    (; β0, β1, β2, β3, β4, b0, b1, b2, b3, b4, b5, b6, b7, a0, a1, a2, a3, a4, a_w_440) = c
    chl = max(chl, oftype(chl, 0.02))           # lower limit of M&M eq. (14)

    # Backscattering, M&M eqs. (11)–(14)
    v = chl > 2 ? zero(chl) : (log10(chl) - b7) / 2
    bbp = (b0 + b1 * (b2 - b3 * log10(chl)) * (λ / b4)^v) * b5 * chl^b6
    bbw = b_w(λ) / 2
    bb = bbw + bbp
    η_b = bbw / bb

    # Absorption, M&M eqs. (16)–(18), with A_chl(440) = 1
    # NOTE: Jin et al.'s supplementary code uses Smith & Baker (1981) a_w below 510 nm but
    # the Pope & Fry (1997) value a_w(440) = 0.00635 in the CDOM term. Kept for parity;
    # Pope & Fry throughout would be consistent with Morel & Maritorena (2001).
    a_y440 = a4 * (a_w_440 + a0 * chl^a1)
    a = a_w(λ) + a0 * A_chl(λ) * chl^a1 + a_y440 * exp(-a2 * (λ - a3))

    # β = β0 - β1η - β2η² + (β3η - β4)μ is linear in μ, so R₀ = A + Bμ
    A = (β0 - β1 * η_b - β2 * η_b^2) * bb / a
    B = (β3 * η_b - β4) * bb / a
    return A, B
end

"""$(TYPEDSIGNATURES) Precompute the coefficients of R₀ = A + Bμ for every band."""
function initialize!(scheme::SpectralOceanVolumeScattering, model::PrimitiveEquation)
    λs = Array(scheme.wavelengths)              # compute on CPU
    A = similar(λs)
    B = similar(λs)
    for k in eachindex(λs)
        A[k], B[k] = R₀_coefficients(scheme, λs[k], scheme.chl)
    end
    copyto!(scheme.R₀_intercept, A)             # works for CPU and GPU arrays
    copyto!(scheme.R₀_slope, B)
    return nothing
end

"""$(TYPEDSIGNATURES) Volume-scattering albedo [1] of band `p.band`, Jin et al. (2011), eq. (6).
In contrast to other `AbstractOceanVolumeScattering`s this needs the band index `band` in `p`, see
`SpectralOceanVolumeScattering` (future work)."""
@inline function (c::SpectralOceanVolumeScattering)(p, μ)
    (; σ, αₛ, band) = p
    (; r0, r1, r2) = c

    r_w = r0 - r1 * σ - r2 * σ^2                                # Jin eq. (7)
    R₀ = c.R₀_intercept[band] + c.R₀_slope[band] * μ            # Jin eqs. (8)–(9), precomputed
    return R₀ * (1 - r_w) * (1 - αₛ) / (1 - r_w * R₀)           # Jin eq. (6)
end