@testset "FFT of geopotential" begin

    # original speedy T30 with 96x48 grid
    P = Parameters(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    geopot_surf_spectral = B.geopot_surf
    geopot_surf_grid = gridded(geopot_surf_spectral,G)

    @test all(geopot_surf_grid .≈ fourier_inverse(fourier(geopot_surf_grid,G),G))
end

@testset "Legendre transform of geopotential" begin

    # original speedy T30 with 96x48 grid
    P = Parameters(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    geopot_surf_spectral = B.geopot_surf
    geopot_surf_grid = gridded(geopot_surf_spectral,G)

    # using the already spectrally truncated geopotential
    @test all(geopot_surf_grid .≈ gridded(spectral(geopot_surf_grid,G),G))
end

@testset "Spectral transform of spectral noise" begin
    # original speedy T30 with 96x48 grid
    P = Parameters(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)

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