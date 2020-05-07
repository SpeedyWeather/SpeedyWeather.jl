using SpeedyWeather
using Test

@testset "FFT of geopotential" begin

    P = Params(T=Float64)
    G = GeoSpectral{P.T}(P)
    B = Boundaries{P.T}(P,G)

    @test all(B.ϕ0trunc .≈ fourier_inverse(fourier(B.ϕ0trunc,G),G))
end

@testset "Legendre transform of geopotential" begin

    P = Params(T=Float64)
    G = GeoSpectral{P.T}(P)
    B = Boundaries{P.T}(P,G)

    # using the already spectrally truncated geopotential
    @test all(B.ϕ0trunc .≈ gridded(spectral(B.ϕ0trunc,G),G))
end

@testset "Spectral transform of spectral noise" begin

    P = Params(T=Float64)
    G = GeoSpectral{P.T}(P)
    B = Boundaries{P.T}(P,G)

    @unpack mx,nx = G.spectral

    A = Complex.(rand(mx,nx))
    At = spectral(gridded(A,G),G)   # the first transform includes truncation

    # the next should be exact to machine precision
    @test all(At .≈ spectral(gridded(At,G),G))
end
