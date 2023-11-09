include("linear_interpolation.jl")

abstract type BandModel end

struct GeneralBandModel{LW, SW, BB, FF} <: BandModel 
    nLW :: Int
    nSW :: Int

    band_limits_LW :: LW
    band_limits_SW :: SW

    blackbody_fraction :: BB
    flux_fraction      :: FF
end

function GeneralBandModel(; limits_SW = ClimLabBandSW(), 
                            limits_LW = optimize_bands(gases, Δν_min = 100.0, κ_min = 1e-6))
    nLW = length(limits_LW)
    nSW = length(limits_SW)

    if limits_SW isa ClimLabBandSW
        flux_fraction = nothing
    else
        flux_fraction = zeros(nSW)
        for band in 1:nSW
            flux_fraction[band] = calc_blackbody_fraction(band, 5777.0, limits_SW)
        end
    end

    if limits_LW isa ClimLabBandLW
        blackbody_fraction = nothing
    else
        temperature = 288-150:0.01:288+80 # Limits of integration to make allow linear interpolation
        blackbody_fraction = []
        for band in 1:nLW
            temp_blackbody = calc_blackbody_fraction.(Ref(band), temperature, Ref(limits_LW))
            push!(blackbody_fraction, LinearInterpolator(temperature, temp_blackbody))
        end
    end

    return GeneralBandModel(nLW, nSW, limits_LW, limits_SW, blackbody_fraction, flux_fraction)
end

bandsizeLW(gbm::GeneralBandModel) = gbm.nLW
bandsizeSW(gbm::GeneralBandModel) = gbm.nSW

blackbodyfraction(band, T, bandmodel) = evaluate(bandmodel.blackbody_fraction[band], T)
fluxfraction(band, bandmodel)         = bandmodel.flux_fraction[band]

calc_blackbody_fraction(band, T, limits) = intensity_fraction(limits[band]..., T)
 
band_limit_SW(gbm::GeneralBandModel) = gbm.band_limits_SW
band_limit_LW(gbm::GeneralBandModel) = gbm.band_limits_LW

function getabsorption_LW(gas, band, bandmodel)
    
    if gas == :O3
        return 0.0
    end

    κ = Symbol(:κ_, gas)
    w = Symbol(:w_, gas)

    limits = band_limit_LW(bandmodel)[band]
    
    @eval begin
        κₙ = absorption_integral($limits..., $w, $κ) 
    end
    return κₙ
end

function getabsorption_SW(gas, band, bandmodel)
    
    κ = Symbol(:κ_, gas)
    w = Symbol(:w_, gas)

    limits = band_limit_SW(bandmodel)[band]
    
    @eval begin
        κₙ = absorption_integral($limits..., $w, $κ) 
    end
    return κₙ
end

function gas_distribution(pressure, gas; T = 288)
    if gas == :O3
        return O3_distribution(pressure)
    elseif gas == :CO2
        return 380e-6 ./ (1 + 380e-6) .* ones(eltype(pressure), size(pressure))
    elseif gas == :H2O
        return H2O_distribution(pressure; T)
    elseif gas == :CH4
        return 1.5e-6 .* ones(eltype(pressure), size(pressure))
    end
end

function H2O_distribution(p; rₛ = 0.77, T = 288)     
    Tcel = @. T - 273.15

    SP = @. 611.2 * exp(17.67 * Tcel / (Tcel + 243.5))
    
    RH = @. rₛ * (p / p₀ - 0.02) / (1 - 0.02)

    VP = RH .* SP

    PP = @. 0.622 * VP / (p - 0.378 * VP)

    PP = max.(PP, Ref(5e-10))

    return PP 
end

function O3_distribution(pressure)

    O3_obs =  [7.82792878e-06, 8.64150529e-06, 7.58940028e-06, 5.24567145e-06,
    3.17761574e-06, 1.82320006e-06, 9.80756960e-07, 6.22870516e-07,
    4.47620550e-07, 3.34481169e-07, 2.62570302e-07, 2.07898125e-07,
    1.57074555e-07, 1.12425545e-07, 8.06004999e-08, 6.27826498e-08,
    5.42990561e-08, 4.99506089e-08, 4.60075681e-08, 4.22977789e-08,
    3.80559071e-08, 3.38768568e-08, 3.12171619e-08, 2.97807119e-08,
    2.87980968e-08, 2.75429934e-08]

    p_obs =   [3.544638,   7.388814,  13.967214,  23.944625,  37.23029 ,  53.114605,
    70.05915  , 85.439115, 100.514695, 118.250335, 139.115395, 163.66207 ,
    192.539935, 226.513265, 266.481155, 313.501265, 368.81798 , 433.895225,
    510.455255, 600.5242  , 696.79629 , 787.70206 , 867.16076 , 929.648875,
    970.55483 , 992.5561  ] .* 100

   return LinearInterpolator(p_obs, O3_obs)(pressure)
end

"""
Simplified band models used by ClimLab
"""
struct ClimLabBandLW end
struct ClimLabBandSW end

length(::ClimLabBandLW) = 4
length(::ClimLabBandSW) = 4

const ClimLabBandLWModel = GeneralBandModel{<:ClimLabBandLW}
const ClimLabBandSWModel = GeneralBandModel{<:Any, <:ClimLabBandSW}

function blackbodyfraction(band, T, ::ClimLabBandLWModel) 
    if band == 1
        return 1.0 - 0.81825 + 3.0E-6 * (T-247.)^2 + 5.225E-6 * (T-282.)^2 - 1.0E-5 *(T-315.)^2
    elseif band == 2
        return 0.148 - 3.0E-6 * (T-247.)^2
    elseif band == 3
        return (0.375 - 5.5E-6 * (T-282.)^2)*0.95
    else
        return 0.314 + 1.0E-5 * (T-315.)^2
    end
end

fluxfraction(band, ::ClimLabBandSWModel) = band == 1 ? 0.06 :
                                           band == 2 ? 0.14 :
                                           band == 3 ? 0.27 : 0.53

absorption_cross_section_LW(::Val{:CO2}, band) = [0.0, 4.0, 0.0, 0.0][band] ./ 1e5 .* g ./ 380e-6
absorption_cross_section_LW(::Val{:H2O}, band) = [0.0, 0.0, 0.7, 50.][band] ./ 1e5 .* g ./ 1e-3
absorption_cross_section_LW(::Val{:O3},  band) = [0.7, 0.0, 0.0, 0.0][band] ./ 1e5 .* g ./ 5e-6

absorption_cross_section_SW(::Val{:CO2}, band) = 0.0
absorption_cross_section_SW(::Val{:H2O}, band) = [0.0, 0.0, 0.0, 0.0012][band]
absorption_cross_section_SW(::Val{:O3},  band) = [200.e-24, 0.0, 0.285e-24, 0.][band] .* gas_constant(Val(:air)) / kB
absorption_cross_section_SW(::Val{:CH4}, band) = 0.0

getabsorption_LW(gas, band, ::ClimLabBandLWModel) = absorption_cross_section_LW(Val(gas), band)
getabsorption_SW(gas, band, ::ClimLabBandSWModel) = absorption_cross_section_SW(Val(gas), band)
