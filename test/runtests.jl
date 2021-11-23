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

@testset "Horizontal diffusion of random" begin
    
    mx,nx,nlev = 31,32,8

    for T in (Float16,Float32,Float64)

        A           = rand(Complex{T},mx,nx,nlev)
        tendency    = rand(Complex{T},mx,nx,nlev)
        D1          = rand(T,mx,nx)
        D2          = rand(T,mx,nx)

        A1 = A[:,:,1]
        tendency1 = copy(tendency[:,:,1])

        SpeedyWeather.do_horizontal_diffusion!(A,tendency,D1,D2)
        SpeedyWeather.do_horizontal_diffusion!(A1,tendency1,D1,D2)

        @test tendency1 == tendency[:,:,1]
    end
end