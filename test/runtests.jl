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
    # original speedy T30 with 96x48 grid
    P = Params(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    mx = G.spectral.mx
    nx = G.spectral.nx

    A = rand(mx,nx) + im*rand(mx,nx)
    At = spectral(gridded(A,G),G)   # the first transform includes truncation

    # the next should be exact to machine precision
    At2 = spectral(gridded(At,G),G)
    At3 = spectral(gridded(At2,G),G)
    @test all(At .≈ At2)
    @test all(At2 .≈ At3)
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

        SpeedyWeather.horizontal_diffusion!(A,tendency,D1,D2)
        SpeedyWeather.horizontal_diffusion!(A1,tendency1,D1,D2)

        @test tendency1 == tendency[:,:,1]
    end
end