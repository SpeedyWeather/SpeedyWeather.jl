@testset "FullGaussianGrid: Test grid and spectral resolution match" begin
    
    p = 5:10

    @testset for dealiasing in [2, 3]
        # powers of two minus 1, T31, T63, T127, etc
        Ts = [2^i - 1 for i in p]                   # spectral resolutions
        nlons = [(dealiasing+1)*2^i for i in p]     # number of longitudes
        nlats = nlons/2                             # number of latitudes
        for (T, nlon, nlat) in zip(Ts, nlons, nlats)
            nlat_half = SpeedyTransforms.get_nlat_half(T, dealiasing)
            @test (nlon, nlat) == (4nlat_half, 2nlat_half)

            trunc = SpeedyTransforms.get_truncation(nlat_half, dealiasing)
            @test T == trunc
        end
    end

    # for these resolutions just test idempotence as the roundup_fft may
    # give various other options than just the 3*2^n-matching
    @testset for dealiasing in [2, 3]
        # T42, T85, T170, T341, T682, T1365 etc
        Ts = [floor(Int, 2^(i+2)/3) for i in p]
        nlons = [dealiasing*2^(i+1) for i in p]
        nlats = nlons/2                             # number of latitudes
        for (T, nlon, nlat) in zip(Ts, nlons, nlats)
            nlat_half = SpeedyTransforms.get_nlat_half(T, dealiasing)
            trunc = SpeedyTransforms.get_truncation(nlat_half, dealiasing)
            nlat_half2 = SpeedyTransforms.get_nlat_half(trunc, dealiasing)
            @test nlat_half == nlat_half2
        end
    end
end

# for the following testsets test some spectral truncations
# but not too large ones as they take so long

spectral_resolutions = (31, 63, 127)
spectral_resolutions_inexact = (127, 255)

@testset "Transform: l=0, m=0 is constant > 0" begin
    for trunc in spectral_resolutions
        for NF in (Float32, Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,
                            HEALPixGrid,
                            OctaHEALPixGrid,
                            FullHEALPixGrid,
                            FullOctaHEALPixGrid)

                SG = SpectralGrid(NF;trunc, Grid)
                S = SpectralTransform(SG)

                alms = zeros(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)
                fill!(alms, 0)
                alms[1, 1] = 1

                map = gridded(alms, S)
            
                for ij in SpeedyWeather.eachgridpoint(map)
                    @test map[ij] ≈ map[1] > zero(NF)
                end
            end
        end
    end
end

@testset "Transform: Recompute, precompute identical results" begin
    for trunc in spectral_resolutions
        for NF in (Float32, Float64)

            SG = SpectralGrid(NF;trunc)
            S1 = SpectralTransform(SG, recompute_legendre=true)
            S2 = SpectralTransform(SG, recompute_legendre=false)

            alms = randn(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)

            map1 = gridded(alms, S1)
            map2 = gridded(alms, S2)
        
            # is only approx as recompute_legendre may use a different precision
            @test map1 ≈ map2
        end
    end
end

@testset "Transform: Individual Legendre polynomials" begin
    @testset for trunc in spectral_resolutions
        for NF in (Float32, Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid)

                SG = SpectralGrid(NF;trunc, Grid)
                S = SpectralTransform(SG, recompute_legendre=true)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)
                        alms[l, m] = 1

                        map = gridded(alms, S)
                        alms2 = spectral(map, S)

                        for lm in SpeedyWeather.eachharmonic(alms, alms2)
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
        @testset for NF in (Float32, Float64)
            @testset for Grid in (  HEALPixGrid,
                                    OctaHEALPixGrid,
                                    FullHEALPixGrid,
                                    FullOctaHEALPixGrid)
                
                SG = SpectralGrid(NF;trunc, Grid)
                S = SpectralTransform(SG, recompute_legendre=true)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)
                        alms[l, m] = 1

                        map = gridded(alms, S)
                        alms2 = spectral(map, S)

                        for lm in SpeedyWeather.eachharmonic(alms, alms2)
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
    @testset for trunc in [31, 42]
        @testset for NF in (Float64, Float32)
            @testset for Grid in (   FullGaussianGrid,
                                     FullClenshawGrid,
                                     OctahedralGaussianGrid,
                                     OctahedralClenshawGrid)

                # clenshaw-curtis grids are only exact for cubic truncation
                dealiasing = Grid in (FullGaussianGrid, OctahedralGaussianGrid) ? 2 : 3

                SG = SpectralGrid(NF;trunc, Grid, dealiasing)
                S = SpectralTransform(SG, recompute_legendre=false)
                O = EarthOrography(SG, smoothing=true, smoothing_truncation=31)
                E = Earth(SG)
                initialize!(O, E, S)

                oro_grid = O.orography
                oro_spec = spectral(oro_grid, S)

                oro_grid1 = gridded(oro_spec, S)
                oro_spec1 = spectral(oro_grid1, S)
                oro_grid2 = gridded(oro_spec1, S)
                oro_spec2 = spectral(oro_grid2, S)

                tol = 1e-1

                for lm in SpeedyWeather.eachharmonic(oro_spec1, oro_spec2)
                    @test oro_spec1[lm] ≈ oro_spec2[lm] atol=tol rtol=tol
                end
                for ij in eachindex(oro_grid1, oro_grid2)
                    @test oro_grid1[ij] ≈ oro_grid2[ij] atol=tol rtol=tol
                end
            end
        end
    end
end
