# The free troposphere layers of the diagnostic clouds are counted on a CPU copy of the σ levels:
# `count` over a view of a concrete array falls back to scalar indexing, which on a device buffer
# (GPU, TPU) is disallowed in non-interactive sessions and would fail `initialize!`

@testset "Free troposphere layers with Reactant vertical coordinates" begin
    clouds = DiagnosticClouds(SpectralGrid(nlayers = 8))
    for VerticalCoordinates in (SigmaCoordinates, SigmaPressureCoordinates), nlayers in (1, 2, 8, 16)
        cpu_coordinates = VerticalCoordinates(SpectralGrid(; nlayers))
        reactant_coordinates = VerticalCoordinates(SpectralGrid(; nlayers, architecture = ReactantDevice()))

        # dispatches to the Reactant extension's CPU-copy method, not the generic one
        method = which(SpeedyWeather.free_troposphere_layers, typeof.((clouds, reactant_coordinates)))
        @test parentmodule(method) == SpeedyWeatherReactantExt

        layers = SpeedyWeather.free_troposphere_layers(clouds, reactant_coordinates)
        @test layers == SpeedyWeather.free_troposphere_layers(clouds, cpu_coordinates)
    end
end
