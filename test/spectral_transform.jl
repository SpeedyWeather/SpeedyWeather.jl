@testset "FFT of geopotential" begin

    # original speedy T30 with 96x48 grid
    P = Params(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    geopot_surf_spectral = B.geopot_surf
    geopot_surf_grid = convert_to_grid(geopot_surf_convert_to_spectral,G)

    @test all(geopot_surf_grid .≈ fourier_inverse(fourier(geopot_surf_grid,G),G))
end

@testset "Legendre transform of geopotential" begin

    # original speedy T30 with 96x48 grid
    P = Params(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)
    B = Boundaries{P.NF}(P,G)

    geopot_surf_spectral = B.geopot_surf
    geopot_surf_grid = convert_to_grid(geopot_surf_spectral,G)

    # using the already spectrally truncated geopotential
    @test all(geopot_surf_grid .≈ convert_to_grid(convert_to_spectral(geopot_surf_grid,G),G))
end

@testset "Spectral transform of spectral noise" begin
    # original speedy T30 with 96x48 grid
    P = Params(NF=Float64;trunc=30,nlat=48,nlon=96)
    G = GeoSpectral{P.NF}(P)

    mx = G.spectral.mx
    nx = G.spectral.nx

    A = rand(mx,nx) + im*rand(mx,nx)
    At = convert_to_spectral(convert_to_grid(A,G),G)   # the first transform includes truncation

    # the next should be exact to machine precision
    At2 = convert_to_spectral(convert_to_grid(At,G),G)
    At3 = convert_to_spectral(convert_to_grid(At2,G),G)
    @test all(At .≈ At2)
    @test all(At2 .≈ At3)
end