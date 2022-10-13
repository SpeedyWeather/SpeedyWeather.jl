@testset "FullGaussianGrid: Test grid and spectral resolution match" begin
    Ts = (31,42,85,170,341,682)     # spectral resolutions
    nlons = (96,128,256,512,1024)   # number of longitudes
    nlats = (48,64,128,256,512)     # number of latitudes
    for (T,nlon,nlat) in zip(Ts,nlons,nlats)
        
        nlat_half = SpeedyWeather.get_resolution(FullGaussianGrid,T)
        @test (nlon,nlat) == (4nlat_half,2nlat_half)

        trunc = SpeedyWeather.get_truncation(FullGaussianGrid,nlat÷2)
        @test T == trunc
    end
end

# for the following testsets test some spectral truncations
# but not too large ones as they take so long
spectral_resolutions = (31,63,127)
spectral_resolutions_inexact = (127,255)

@testset "Transform: l=0,m=0 is constant > 0" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,
                            HEALPixGrid,
                            HEALPix4Grid,
                            FullHEALPixGrid,
                            FullHEALPix4Grid)

                p,d,m = initialize_speedy(NF;trunc,Grid)
                S = m.spectral_transform

                alms = copy(p.layers[1].leapfrog[1].vor)
                fill!(alms,0)
                alms[1,1] = 1

                map = gridded(alms,S)
            
                for ij in SpeedyWeather.eachgridpoint(map)
                    @test map[ij] ≈ map[1] > zero(NF)
                end
            end
        end
    end
end

@testset "Transform: Recompute, precompute identical results" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            p1,d1,m1 = initialize_speedy(NF;trunc,recompute_legendre=false)
            p2,d2,m2 = initialize_speedy(NF;trunc,recompute_legendre=true)

            (;vor) = p1.layers[1].leapfrog[1]
            alms = randn(typeof(vor),size(vor)...)

            map1 = gridded(alms,m1.spectral_transform)
            map2 = gridded(alms,m2.spectral_transform)
        
            # is only approx as recompute_legendre may use a different precision
            @test map1 ≈ map2
        end
    end
end

@testset "Transform: Individual Legendre polynomials" begin
    @testset for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid)

                P = Parameters(;NF,trunc,Grid)
                S = SpectralTransform(P)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}},S.lmax+2,S.mmax+1)
                        alms[l,m] = 1

                        map = gridded(alms,S)
                        alms2 = spectral(map,S)

                        for lm in SpeedyWeather.eachharmonic(alms,alms2)
                            @test alms[lm] ≈ alms2[lm] atol=100*eps(NF)
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Individual Legendre polynomials (inexact transforms)" begin
    @testset for trunc in spectral_resolutions_inexact
        @testset for NF in (Float32,Float64)
            @testset for Grid in (  HEALPixGrid,
                                    HEALPix4Grid,
                                    FullHEALPixGrid,
                                    FullHEALPix4Grid)
                P = Parameters(;NF,trunc,Grid)
                S = SpectralTransform(P)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}},S.lmax+2,S.mmax+1)
                        alms[l,m] = 1

                        map = gridded(alms,S)
                        alms2 = spectral(map,S)

                        for lm in SpeedyWeather.eachharmonic(alms,alms2)
                            @test alms[lm] ≈ alms2[lm] atol=1e-3 rtol=1e-3
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Orography (exact grids)" begin

    # Test for variable resolution
    @testset for trunc in [31,42]
        @testset for NF in (Float64,Float32)
            @testset for Grid in (   FullGaussianGrid,
                                     FullClenshawGrid,
                                     OctahedralGaussianGrid,
                                     OctahedralClenshawGrid)
                            
                P = Parameters(;NF,Grid,trunc,model=:shallowwater)
                S = SpectralTransform(P)
                B = Boundaries(P,S)

                oro_grid = B.orography
                oro_spec = spectral(oro_grid,S)

                # smooth orography
                for m in 1:trunc+1
                    for l in m:trunc+2
                        oro_spec[l,m] = 0
                    end
                end 

                oro_grid1 = gridded(oro_spec,S)
                oro_spec1 = spectral(oro_grid1,S)
                oro_grid2 = gridded(oro_spec1,S)
                oro_spec2 = spectral(oro_grid2,S)

                tol = 1e-1

                for lm in SpeedyWeather.eachharmonic(oro_spec1,oro_spec2)
                    @test oro_spec1[lm] ≈ oro_spec2[lm] atol=tol rtol=tol
                end
                for ij in eachindex(oro_grid1,oro_grid2)
                    @test oro_grid1[ij] ≈ oro_grid2[ij] atol=tol rtol=tol
                end
            end
        end
    end
end
