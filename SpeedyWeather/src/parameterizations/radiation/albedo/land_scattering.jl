abstract type AbstractLandScattering <: AbstractScattering end

"""Visible radiation isotropic scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Isotropic scattering means that
light at the surface is reflected equally in all directions, regardless of where the sun is.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleIsotropicScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.13926611446264833

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.049407226600965096

    "Coefficient fit via symbolic regression"
    @param c3::NF = 0.03012000138970739

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.007003960811839407

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.0030819533608978183

    "Coefficient fit via symbolic regression"
    @param c6::NF = 98.73431990602951

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.0007767787214513013

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.4746273568956652

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.004805110395235223

    "Coefficient fit via symbolic regression"
    @param c10::NF = 1.1381888960412316

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.0203040824498551

end

Adapt.@adapt_structure VisibleIsotropicScattering
VisibleIsotropicScattering(SG::SpectralGrid; kwargs...) = VisibleIsotropicScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleIsotropicScattering)(p)
    (; laiₗ, sm1, sm2, snow_cover, st2, t2m, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = c

    # Mean calculation
    μ = (
        -(c1 * vₕ) + vₗ * (-(c2 * laiₗ) - c2 * snow_cover - c2 * sm1 - c3) + snow_cover * (c4 - c5 * vₕ)
            * (ΔTₛ + c6) - c1 * sm2 - c7 * t2m + c8
    )

    # Std deviation
    raw_σ = snow_cover * sm2 * (c9 * st2 - c10) + c11
    σ = max(raw_σ, zero(raw_σ))  # shouldn't be negative...
    return μ + σ * noise
end

"""Visible radiation volumetric scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Volumetric scattering represents the
reflection mechanism of light interacting with horizontally homogenous vegetation canopy at the surface - 
light meets the canopy, bounces around a bit beneath the canopy and exits back towards the sky. This 
does depend on where the sun is and the viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleVolumetricScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.005137089001311701

    "Coefficient fit via symbolic regression"
    @param c2::NF = 8.359579016151135e-7

    "Coefficient fit via symbolic regression"
    @param c3::NF = 547.1867710914254

    "Coefficient fit via symbolic regression"
    @param c4::NF = 4.88214507804199

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.006321233141025383

    "Coefficient fit via symbolic regression"
    @param c6::NF = 547.1858817550591

    "Coefficient fit via symbolic regression"
    @param c7::NF = 7.789539669853218e-8

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.0012984891007226625

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.013366006624365259

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.037107779858204085

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.001076347986980297

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.02154070145308177

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.013735265859374146

end

Adapt.@adapt_structure VisibleVolumetricScattering
VisibleVolumetricScattering(SG::SpectralGrid; kwargs...) = VisibleVolumetricScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleVolumetricScattering)(p)
    (; geopotential, laiₗ, sm1, snow_cover, st1, t2m, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c
    # Vegetation cover term
    q0 = vₕ + vₗ

    # Mean calculation
    μ = (
        q0 * (
            -(c1 * laiₗ) + c2 * t2m * (
                ΔTₛ + snow_cover * (
                    -(c3 * snow_cover) + c4 * t2m + c5 *
                        geopotential - c6
                ) + c7 * geopotential^2
            ) / (q0 + c8) - c9
        ) + c8 * st1 - c10 * sm1 - c11 * t2m
    )

    # Std devation calculation
    raw_σ = c12 * vₗ * snow_cover + c13
    σ = max(raw_σ, zero(raw_σ))  # can't be negative
    return μ + σ * noise
end

"""Visible radiation geometric-optical scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Geometric scattering represents the 
reflection due to shadowcasting of 3D objects at the surface, like tree crowns or big rocks, and 
obscuration of shadows due to viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleGeometricScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 98.30689271056275

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.0002465039456962893

    "Coefficient fit via symbolic regression"
    @param c3::NF = 1.74846262394399e-8

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.2655102008045942

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.00016722961466941218

    "Coefficient fit via symbolic regression"
    @param c6::NF = 1.8198639960846736e-6

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.009976526656995078

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.000145787380943221

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.012478157859566873

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.01037571803805248

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.007018418642708496

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.000694635309383557

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.005474843497829615

end

Adapt.@adapt_structure VisibleGeometricScattering
VisibleGeometricScattering(SG::SpectralGrid; kwargs...) = VisibleGeometricScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleGeometricScattering)(p)
    (; U₁₀, geopotential, laiₗ, sm1, snow_cover, snow_depth, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c

    # Mean calculation
    μ = (
        snow_cover * (snow_depth * (ΔTₛ - U₁₀) + c1) * (
            -(c2 * vₗ) + c3 * geopotential * (
                vₕ + laiₗ * (vₗ + snow_depth) - c4
            ) + c5
        ) + (
            c6 * geopotential + (-(c7 * vₕ) - c7 * vₗ) * (
                laiₗ * sm1 + c8 * geopotential
            ) + c9
        ) * sqrt(vₗ) + c10
    )

    # Std deviation
    raw_σ = (
        snow_cover * (c11 - c12 * snow_depth) + c13
    )
    σ = max(raw_σ, zero(raw_σ))  # can't be negative
    return μ + σ * noise
end

"""Near-infrared radiation isotropic scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Isotropic scattering means that
light at the surface is reflected equally in all directions, regardless of where the sun is.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRIsotropicScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.13250859400402135

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.08409939463271902

    "Coefficient fit via symbolic regression"
    @param c3::NF = 0.03616646161624834

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.10051756565160144

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.200655651610176

    "Coefficient fit via symbolic regression"
    @param c6::NF = 0.13125739630648842

    "Coefficient fit via symbolic regression"
    @param c7::NF = 1.2330577987093237e-5

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.308490327357136

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.38361527926198996

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.5104495709936557

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.07677930381946513

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.052701989554117486

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.25930825098307064

    "Coefficient fit via symbolic regression"
    @param c14::NF = 0.047761528867464

    "Coefficient fit via symbolic regression"
    @param c15::NF = 0.29940528501431796

    "Coefficient fit via symbolic regression"
    @param c16::NF = 0.03240815843538468

end

Adapt.@adapt_structure NIRIsotropicScattering
NIRIsotropicScattering(SG::SpectralGrid; kwargs...) = NIRIsotropicScattering{SG.NF}(; kwargs...)

@inline function (c::NIRIsotropicScattering)(p)
    (; geopotential, laiₗ, sm2, snow_cover, vₕ, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16) = c

    # Mean calculation
    μ = (
        vₕ * (-(c1 * snow_cover) - c2) + laiₗ * (-(c3 * vₗ) + c4 * sm2) + c5 * snow_cover^3 + (
            - (c6 * laiₗ) + c7 * geopotential * (c8 * snow_cover + sm2 - c9)
        ) / (laiₗ + snow_cover + c10) - (c11 * sm2 - c12) / sqrt(c13 * sm2 + c14) + c15
    )

    # Std deviation, in this case a constant
    σ = c16
    return μ + σ * noise
end

"""Near-infrared radiation volumetric scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Volumetric scattering represents the
reflection mechanism of light interacting with horizontally homogenous vegetation canopy at the surface - 
light meets the canopy, bounces around a bit beneath the canopy and exits back towards the sky. This 
does depend on where the sun is and the viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRVolumetricScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.02644856368910439

    "Coefficient fit via symbolic regression"
    @param c2::NF = 2.7037305935790665

    "Coefficient fit via symbolic regression"
    @param c3::NF = 2.2973971005827364e-5

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.0063074532998523255

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.015870141431257815

    "Coefficient fit via symbolic regression"
    @param c6::NF = 7.900349778266809e-8

    "Coefficient fit via symbolic regression"
    @param c7::NF = 1.3580861013567708e-11

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.1966651392522715

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.026503362065879634

end

Adapt.@adapt_structure NIRVolumetricScattering
NIRVolumetricScattering(SG::SpectralGrid; kwargs...) = NIRVolumetricScattering{SG.NF}(; kwargs...)

@inline function (c::NIRVolumetricScattering)(p)
    (; sm2, snow_cover, st1, st2, t2m, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9) = c

    # Mean calculation
    μ = (
        snow_cover * (c1 * vₗ + c2) + sm2 * (c3 * st1^2 - c4 * st2) + t2m * (
            - (c5 * snow_cover) + t2m^2 * (c6 * snow_cover - c7 * st1)
        ) + c8
    )

    # Std deviation, in this case just a constant
    σ = c9
    return μ + σ * noise
end

"""Near-infrared radiation geometric-optical scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Geometric scattering represents the 
reflection due to shadowcasting of 3D objects at the surface, like tree crowns or big rocks, and 
obscuration of shadows due to viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRGeometricScattering{NF} <: AbstractLandScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.00016970668634889394

    "Coefficient fit via symbolic regression"
    @param c2::NF = 4.1657057100227406e-10

    "Coefficient fit via symbolic regression"
    @param c3::NF = 1.8792166707362004e-7

    "Coefficient fit via symbolic regression"
    @param c4::NF = 2.0653707823044138e-5

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.04338449131472477

    "Coefficient fit via symbolic regression"
    @param c6::NF = 0.10172346530005257

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.07360077673775481

    "Coefficient fit via symbolic regression"
    @param c8::NF = 3.49400847746966e-5

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.005699938844269787

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.0008285983510980352

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.614697095826219

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.0790756374361683

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.009887490927614274

end

Adapt.@adapt_structure NIRGeometricScattering
NIRGeometricScattering(SG::SpectralGrid; kwargs...) = NIRGeometricScattering{SG.NF}(; kwargs...)

@inline function (c::NIRGeometricScattering)(p)
    (; geopotential, laiₕ, sm2, snow_cover, st2, t2m, vₕ, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c

    # Mean calculation
    μ = (
        -(c1 * t2m) + geopotential * (c2 * st2^2 - c3 * st2 + c4) + (c5 - c6 * snow_cover) * (
            vₕ * (
                c7 * laiₕ - c8 * geopotential
            ) + sm2 - c9 * t2m * (-(c10 * st2) + sm2) * (vₗ + 2 * sm2 - c11)
        ) + c12
    )

    # Std deviation, in this case a constant
    σ = c13
    return μ + σ * noise
end