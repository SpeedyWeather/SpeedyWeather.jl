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
                            OctaminimalGaussianGrid,
                            HEALPixGrid,
                            OctaHEALPixGrid,
                            FullHEALPixGrid,
                            FullOctaHEALPixGrid)

                SG = SpectralGrid(; NF, trunc, Grid)
                S = SpectralTransform(SG)

                alms = zeros(LowerTriangularMatrix{Complex{NF}}, SG.spectrum)
                fill!(alms, 0)
                alms[1, 1] = 1

                map = transform(alms, S)
            
                for ij in SpeedyWeather.eachgridpoint(map)
                    @test map[ij] ≈ map[1] > zero(NF)
                end
            end
        end
    end
end

# functionality deprecated
# @testset "Transform: Recompute, precompute identical results" begin
#     for trunc in spectral_resolutions
#         for NF in (Float32, Float64)

#             SG = SpectralGrid(; NF, trunc)
#             S1 = SpectralTransform(SG, recompute_legendre=true)
#             S2 = SpectralTransform(SG, recompute_legendre=false)

#             alms = randn(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)

#             map1 = transform(alms, S1)
#             map2 = transform(alms, S2)
        
#             # is only approx as recompute_legendre may use a different precision
#             @test map1 ≈ map2
#         end
#     end
# end

@testset "Transform: Individual Legendre polynomials" begin
    @testset for trunc in spectral_resolutions
        @testset for NF in (Float32, Float64)
            @testset for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid)

                SG = SpectralGrid(; NF, trunc, Grid)
                S = SpectralTransform(SG)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, SG.spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)
                        alms2 = transform(map, S)

                        for lm in SpeedyWeather.eachharmonic(alms, alms2)
                            @test alms[lm] ≈ alms2[lm] atol=100*eps(NF)
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Singleton dimensions" begin
    @testset for trunc in spectral_resolutions
        for NF in (Float32, Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,
                            OctaminimalGaussianGrid)

                SG = SpectralGrid(; NF, trunc, Grid)
                S = SpectralTransform(SG)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, SG.spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)
                        
                        # add singleton dimension for lower triangular matrix
                        alms2 = zeros(LowerTriangularArray{Complex{NF}}, SG.spectrum, 1)
                        alms2[:, 1] = alms
                        map2 = deepcopy(map)
                        transform!(map2, alms2, S)

                        for ij in eachindex(map, map2)
                            @test map[ij] == map2[ij]
                        end

                        # add singleton dimension for grid
                        grid = randn(SG.GridVariable2D, SG.nlat_half)
                        alms = transform(grid, S)
                        alms2 = deepcopy(alms)

                        grid3D = zeros(SG.GridVariable3D, SG.nlat_half, 1)
                        grid3D[:, 1] = grid

                        transform!(alms2, grid3D, S)

                        for lm in eachindex(alms, alms2)
                            @test alms[lm] == alms2[lm]
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Real to real transform" begin
    for NF in (Float32, Float64)
        # test Float64 -> Float32
        spectral_grid = SpectralGrid(; NF)

        # create real and complex otherwise identical
        Lreal = randn(LowerTriangularMatrix{NF}, spectral_grid.spectrum)
        Lcomplex = complex.(Lreal)
        
        grid1 = zeros(spectral_grid.GridVariable2D, spectral_grid.nlat_half)
        grid2 = zeros(spectral_grid.GridVariable2D, spectral_grid.nlat_half)
        
        S = SpectralTransform(spectral_grid)

        transform!(grid1, Lreal, S)
        transform!(grid2, Lcomplex, S)

        for ij in eachindex(grid1, grid2)
            @test grid1[ij] ≈ grid2[ij]
        end
    end
end

@testset "Transform: NF flexibility spectral inputs" begin

    # test Float64 -> Float32
    spectral_grid = SpectralGrid()
    L_f32 = randn(LowerTriangularMatrix{ComplexF32}, spectral_grid.spectrum)
    L_f64 = ComplexF64.(L_f32)
    
    grid1 = zeros(spectral_grid.GridVariable2D, spectral_grid.nlat_half)
    grid2 = zeros(spectral_grid.GridVariable2D, spectral_grid.nlat_half)
    
    S = SpectralTransform(spectral_grid)

    transform!(grid1, L_f64, S)  # Float64 -> Float32
    transform!(grid2, L_f32, S)  # Float32 -> Float32

    for ij in eachindex(grid1, grid2)
        @test grid1[ij] ≈ grid2[ij] atol=sqrt(eps(Float32))
    end

    # TODO NF flexibility for grid inputs currently not supported for a good reason:
    # the fourier transform is pre-planned for number type NF so you can't use another one
    # above in spectral->grid this doesn't matter as the fourier transform is applied after
    # the Legendre transform (which isn't pre-planned and hence NF flexible)
    # consequently, the following is commented out

    # L1 = deepcopy(L_f32)
    # L2 = deepcopy(L_f32)
    # fill!(L1, 0)
    # fill!(L2, 0)

    # grid3_f32 = randn(spectral_grid.GridVariable2D, spectral_grid.nlat_half)
    # grid3_f64 = Float64.(grid3_f32)

    # transform!(L1, grid3_f32, S)
    # transform!(L2, grid3_f64, S)    # throws an error
end

@testset "Transform: NF flexibility for spectral outputs" begin

    # test Float32 -> Float64
    spectral_grid = SpectralGrid()
    NF = spectral_grid.NF
    grid = randn(spectral_grid.Grid{NF}, spectral_grid.nlat_half)
    
    L1 = zeros(LowerTriangularMatrix{Complex{Float32}}, spectral_grid.spectrum)
    L2 = zeros(LowerTriangularMatrix{Complex{Float64}}, spectral_grid.spectrum)
    
    S = SpectralTransform(spectral_grid)

    transform!(L1, grid, S)  # Float32 -> Float32
    transform!(L2, grid, S)  # Float32 -> Float64

    for lm in eachindex(L1, L2)
        @test L1[lm] ≈ L2[lm]
    end

    # TODO similar to the above, NF flexibility for grid outputs currently not supported because
    # the fourier transform is pre-planned and inplace with LinearAlgebra.mul! which cannot write into
    # another eltype than the FFT plan
end

@testset "Transform: Individual Legendre polynomials (inexact transforms)" begin
    @testset for trunc in spectral_resolutions_inexact
        @testset for NF in (Float32, Float64)
            @testset for Grid in (  HEALPixGrid,
                                    OctaHEALPixGrid,
                                    OctaminimalGaussianGrid,
                                    FullHEALPixGrid,
                                    FullOctaHEALPixGrid)
                
                SG = SpectralGrid(; NF, trunc, Grid)
                S = SpectralTransform(SG,)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, SG.spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)
                        alms2 = transform(map, S)

                        tol = 1e-3

                        for lm in SpeedyWeather.eachharmonic(alms, alms2)
                            @test alms[lm] ≈ alms2[lm] atol=tol rtol=tol
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

                SG = SpectralGrid(; NF, trunc, Grid, dealiasing)
                S = SpectralTransform(SG)
                O = EarthOrography(SG, smoothing=true)
                E = Earth(SG)
                initialize!(O, E, S)

                oro_grid = O.orography
                oro_spec = transform(oro_grid, S)

                oro_grid1 = transform(oro_spec, S)
                oro_spec1 = transform(oro_grid1, S)
                oro_grid2 = transform(oro_spec1, S)
                oro_spec2 = transform(oro_grid2, S)

                tol = 1e-3

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
