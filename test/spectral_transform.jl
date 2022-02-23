@testset "FFT of geopotential" begin

    # Test for variable resolution
    for trunc in (30,42,63,85,127,255,511)
    
        P = Parameters(NF=Float64;trunc=trunc)
        G = GeoSpectral{P.NF}(P)
        B = Boundaries{P.NF}(P,G)

        geopot_surf_spectral = B.geopot_surf
        geopot_surf_grid = gridded(geopot_surf_spectral,G)

        @test all(geopot_surf_grid .≈ fourier_inverse(fourier(geopot_surf_grid,G),G))
    end
end

@testset "Legendre transform of geopotential" begin

    # Test for variable resolution
    for trunc in (30,42,63,85,127,255,511)
    
        P = Parameters(NF=Float64;trunc=trunc)
        G = GeoSpectral{P.NF}(P)
        B = Boundaries{P.NF}(P,G)

        geopot_surf_spectral = B.geopot_surf
        geopot_surf_grid = gridded(geopot_surf_spectral,G)

        # using the already spectrally truncated geopotential
        @test geopot_surf_grid ≈ gridded(spectral(geopot_surf_grid,G),G)
    end
end

@testset "Spectral transform of spectral noise" begin
    
    # Test for variable resolution
    for trunc in (30,42,63,85,127,255,511)
    
        P = Parameters(NF=Float64;trunc=trunc)
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
end