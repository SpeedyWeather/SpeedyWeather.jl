# Output variables used to be defined as e.g. `SpeedyWeather.TemperatureOutput`, they are now
# `OutputVariables.Temperature`. Forward the old names with a deprecation warning
const DEPRECATED_OUTPUT_VARIABLES = (
    :AlbedoOutput => :Albedo,
    :AllOutputVariables => :All,
    :BoundaryLayerDragOutput => :BoundaryLayerDrag,
    :BoundaryLayerOutput => :BoundaryLayer,
    :BoundaryOutput => :Boundary,
    :CloudCoverOutput => :CloudCover,
    :CloudTopOutput => :CloudTop,
    :ConvectiveRainOutput => :ConvectiveRain,
    :ConvectiveRainRateOutput => :ConvectiveRainRate,
    :ConvectiveSnowOutput => :ConvectiveSnow,
    :ConvectiveSnowRateOutput => :ConvectiveSnowRate,
    :DivergenceOutput => :Divergence,
    :DynamicsOutput => :Dynamics,
    :HumidityOutput => :Humidity,
    :InterfaceDisplacementOutput => :InterfaceDisplacement,
    :LandOutput => :Land,
    :LandSeaMaskOutput => :LandSeaMask,
    :LandSurfaceRoughnessOutput => :LandSurfaceRoughness,
    :LargeScaleRainOutput => :LargeScaleRain,
    :LargeScaleRainRateOutput => :LargeScaleRainRate,
    :LargeScaleSnowOutput => :LargeScaleSnow,
    :LargeScaleSnowRateOutput => :LargeScaleSnowRate,
    :MeanSeaLevelPressureOutput => :MeanSeaLevelPressure,
    :MeridionalVelocity10mOutput => :MeridionalVelocity10m,
    :MeridionalVelocityOutput => :MeridionalVelocity,
    :OceanOutput => :Ocean,
    :OceanSurfaceRoughnessOutput => :OceanSurfaceRoughness,
    :OrographyOutput => :Orography,
    :OutgoingLongwaveRadiationOutput => :OutgoingLongwaveRadiation,
    :OutgoingShortwaveRadiationOutput => :OutgoingShortwaveRadiation,
    :PrecipitationOutput => :Precipitation,
    :RadiationOutput => :Radiation,
    :RandomPatternOutput => :RandomPattern,
    :RiverRunoffOutput => :RiverRunoff,
    :SeaIceConcentrationOutput => :SeaIceConcentration,
    :SeaSurfaceTemperatureOutput => :SeaSurfaceTemperature,
    :SnowDepthOutput => :SnowDepth,
    :SnowMeltOutput => :SnowMelt,
    :SoilMoistureAvailabilityOutput => :SoilMoistureAvailability,
    :SoilMoistureOutput => :SoilMoisture,
    :SoilTemperatureOutput => :SoilTemperature,
    :SurfaceFluxesOutput => :SurfaceFluxes,
    :SurfaceHumidityFluxOutput => :SurfaceHumidityFlux,
    :SurfaceLatentHeatFluxOutput => :SurfaceLatentHeatFlux,
    :SurfaceLongwaveDownOutput => :SurfaceLongwaveDown,
    :SurfaceLongwaveUpOutput => :SurfaceLongwaveUp,
    :SurfacePressureOutput => :SurfacePressure,
    :SurfaceRoughnessOutput => :SurfaceRoughness,
    :SurfaceSensibleHeatFluxOutput => :SurfaceSensibleHeatFlux,
    :SurfaceShortwaveDownOutput => :SurfaceShortwaveDown,
    :SurfaceShortwaveUpOutput => :SurfaceShortwaveUp,
    :SurfaceTemperatureOutput => :SurfaceTemperature,
    :TemperatureOutput => :Temperature,
    :TerrariumOutput => :Terrarium,
    :TracerOutput => :Tracer,
    :VorticityOutput => :Vorticity,
    :ZonalVelocity10mOutput => :ZonalVelocity10m,
    :ZonalVelocityOutput => :ZonalVelocity,
)

for (old, new) in DEPRECATED_OUTPUT_VARIABLES
    docstring = "Deprecated, use [`OutputVariables.$new`](@ref) instead."
    message = "`SpeedyWeather.$old` is deprecated, use `OutputVariables.$new` instead."
    @eval begin
        @doc $docstring
        function $old(args...; kwargs...)
            Base.depwarn($message, $(QuoteNode(old)); force = true)
            return OutputVariables.$new(args...; kwargs...)
        end
    end
end
