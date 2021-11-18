using SpeedyWeather
using Test

@testset "FFT of geopotential" begin

    P = Params(NF=Float64)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    @test all(B.ϕ0trunc .≈ fourier_inverse(fourier(B.ϕ0trunc,G),G))
end

@testset "Legendre transform of geopotential" begin

    P = Params(NF=Float64)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    # using the already spectrally truncated geopotential
    @test all(B.ϕ0trunc .≈ gridded(spectral(B.ϕ0trunc,G),G))
end

@testset "Spectral transform of spectral noise" begin

    P = Params(NF=Float64)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    mx = G.spectral.mx
    nx = G.spectral.nx

    A = Complex.(rand(mx,nx))
    At = spectral(gridded(A,G),G)   # the first transform includes truncation

    # the next should be exact to machine precision
    @test all(At .≈ spectral(gridded(At,G),G))
end
