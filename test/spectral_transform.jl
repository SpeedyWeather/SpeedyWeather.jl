# using Test
# include("../src/SpeedyWeather.jl")
# using .SpeedyWeather

spectral_resolutions = (31,42,85,170,341)

@testset "Transform: l=0,m=0 is constant > 0" begin

    # Test for variable resolution
    for trunc in spectral_resolutions
    
        NF = Float64
        P = Parameters(recompute_legendre=true;NF,trunc)
        G = GeoSpectral(P)
        S = G.spectral

        alms = zeros(Complex{NF},S.lmax+1,S.mmax+1)
        alms[1,1] = 1

        map = gridded(alms,S)
        
        for i in eachindex(map)
            @test map[i] == map[1] > zero(NF)
        end
    end
end

@testset "Transform: Recompute, precompute identical results" begin

    # Test for variable resolution
    for trunc in spectral_resolutions
    
        NF = Float64
        P1 = Parameters(recompute_legendre=true;NF,trunc)
        G1 = GeoSpectral(P1)
        S1 = G1.spectral

        P2 = Parameters(recompute_legendre=false;NF,trunc)
        G2 = GeoSpectral(P2)
        S2 = G2.spectral

        alms = rand(Complex{NF},S1.lmax+1,S1.mmax+1)
        map1 = gridded(alms,S1)
        map2 = gridded(alms,S2)
        
        @test map1 == map2
    end
end

@testset "Transform: Individual Legendre polynomials" begin

    # Test for variable resolution
    for trunc in spectral_resolutions
    
        NF = Float64
        P = Parameters(;NF,trunc)
        G = GeoSpectral(P)
        S = G.spectral

        lmax = 10
        for l in 1:lmax
            for m in 1:l
                alms = zeros(Complex{NF},S.lmax+1,S.mmax+1)
                alms[l,m] = 1

                map = gridded(alms,S)
                alms2 = spectral(map,S)

                for i in 1:S.lmax
                    for j in 1:i
                        @test alms[i,j] ≈ alms2[i,j] atol=1e-3 skip=true
                    end
                end
            end
        end
    end
end

@testset "Transform: Geopotential" begin

    # Test for variable resolution
    for trunc in spectral_resolutions
    
        P = Parameters(NF=Float64;trunc=trunc)
        G = GeoSpectral(P)
        B = Boundaries(P,G)
        S = G.spectral

        geopot_surf_spectral = B.geopot_surf
        geopot_surf_grid = gridded(geopot_surf_spectral,S)
        geopot_surf_spectral2 = spectral(geopot_surf_grid,S)
        geopot_surf_grid2 = gridded(geopot_surf_spectral2,S)

        @test all(geopot_surf_spectral .≈ geopot_surf_spectral2) broken=true
        @test all(geopot_surf_grid .≈ geopot_surf_grid2) broken=true
    end
end