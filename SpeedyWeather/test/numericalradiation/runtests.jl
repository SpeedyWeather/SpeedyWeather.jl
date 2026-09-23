# Basic functionality of the NumericalRadiation extension. The full unit tests and the
# validation of the coupling live in NumericalRadiation.jl (test/speedyweather/, validation/).
using SpeedyWeather
using NumericalRadiation
using Test

@testset "NumericalRadiation extension" begin
    include("basic.jl")
end
