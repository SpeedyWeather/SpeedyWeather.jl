module SpeedyWeatherNumericalRadiationExt

# Radiation schemes of NumericalRadiation.jl as SpeedyWeather components. The scheme types are
# NumericalRadiation's own (`ClearSkyEcCKDRadiation`, `AnalyticBandLongwave`); this extension adds the
# methods SpeedyWeather calls (`variables`, `initialize!`, `parameterization!`) and constructors
# from a SpectralGrid. Both are used as any other scheme, e.g.
#   PrimitiveWetModel(spectral_grid; radiation = ClearSkyEcCKDRadiation(spectral_grid))
#   Radiation(spectral_grid; longwave = AnalyticBandLongwave(spectral_grid))

using SpeedyWeather
using NumericalRadiation
using DocStringExtensions

import NumericalRadiation: AtmosphereProfile, ColumnGrid, SurfaceState,
    PhysicalConstants, LongwaveDiagnostics, solve_longwave!, AnalyticBandLongwave
import NumericalRadiation: ClearSkyEcCKDRadiation, EcCKDTabulatedGasOpticsModel, ColumnAtmosphere, RadiativeFluxes,
    LongwaveOptics, ShortwaveOptics, CloudlessLongwave, CloudlessShortwave,
    ShortwaveColumnScratch, TabulatedSurfaceEmission, LongwaveBoundaryConditions,
    ShortwaveBoundaryConditions, optical_properties!, radiative_fluxes!,
    read_reference_ecckd_gas_optics

# CO₂ [ppm] of the analytic-band longwave when the model has no `greenhouse_gases.co2`;
# NumericalRadiation's own default (AtmosphereProfile). ClearSkyEcCKDRadiation carries its own.
const DEFAULT_CO₂ = 280

# The schemes are NumericalRadiation's types, not SpeedyWeather model components. The only
# place SpeedyWeather dispatches on the component type is the time-step selection
# (get_step, get_prognostic_step, get_tendency_step); for that the schemes take the steps
# any radiation parameterization gets, forwarded through a stand-in of that type.
struct RadiationStandIn <: SpeedyWeather.AbstractRadiation end
const NumericalRadiationScheme = Union{AnalyticBandLongwave, ClearSkyEcCKDRadiation}
for f in (:get_step, :get_prognostic_step, :get_tendency_step)
    @eval begin
        @inline SpeedyWeather.$f(var, TS::SpeedyWeather.AbstractTimeStepper, ::NumericalRadiationScheme) =
            SpeedyWeather.$f(var, TS, RadiationStandIn())
        @inline SpeedyWeather.$f(var, TS::SpeedyWeather.AbstractTimeStepper, ::NumericalRadiationScheme, M::SpeedyWeather.AbstractModel) =
            SpeedyWeather.$f(var, TS, RadiationStandIn(), M)
    end
end

include("analytic_band_longwave.jl")   # AnalyticBandLongwave as Radiation(; longwave)
include("ecckd_radiation.jl")          # ClearSkyEcCKDRadiation: clear-sky ecCKD, both streams in one component

end # module
