@testset "FullGaussianGrid: Test grid and spectral resolution match" begin

    p = 5:10

    @testset for dealiasing in [2, 3]
        # powers of two minus 1, T31, T63, T127, etc
        Ts = [2^i - 1 for i in p]                   # spectral resolutions
        nlons = [(dealiasing + 1) * 2^i for i in p]     # number of longitudes
        nlats = nlons / 2                             # number of latitudes
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
        Ts = [floor(Int, 2^(i + 2) / 3) for i in p]
        nlons = [dealiasing * 2^(i + 1) for i in p]
        nlats = nlons / 2                             # number of latitudes
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
    @testset for trunc in spectral_resolutions
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    FullClenshawGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid,
                )

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                S = SpectralTransform(spectrum, grid; NF)

                alms = zeros(LowerTriangularMatrix{Complex{NF}}, spectrum)
                fill!(alms, 0)
                alms[1, 1] = 1

                map = transform(alms, S)

                for ij in eachgridpoint(map)
                    @test map[ij] ≈ map[1] > zero(NF)
                end
            end
        end
    end
end

@testset "Transform: Individual Legendre polynomials" begin
    @testset for trunc in spectral_resolutions
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    FullClenshawGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                )

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                S = SpectralTransform(spectrum, grid; NF)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)
                        alms2 = transform(map, S)

                        for lm in eachharmonic(alms, alms2)
                            @test alms[lm] ≈ alms2[lm] atol = 100 * eps(NF)
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Singleton dimensions" begin
    @testset for trunc in spectral_resolutions
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    FullClenshawGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
                )

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                S = SpectralTransform(spectrum, grid; NF)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}}, spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)

                        # add singleton dimension for lower triangular matrix
                        alms2 = zeros(Complex{NF}, spectrum, 1)
                        alms2[:, 1] = alms
                        map2 = deepcopy(map)
                        transform!(map2, alms2, S)

                        for ij in eachindex(map, map2)
                            @test map[ij] == map2[ij]
                        end

                        # add singleton dimension for grid
                        field = randn(NF, grid)
                        alms = transform(field, S)
                        alms2 = deepcopy(alms)

                        field3D = zeros(NF, grid, 1)
                        field3D[:, 1] = field       # without broadcasting
                        field3D[:, 1] .= field      # with

                        transform!(alms2, field3D, S)

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
        spectrum = Spectrum(31)
        grid = OctahedralGaussianGrid(24)
        S = SpectralTransform(spectrum, grid; NF)

        # create real and complex otherwise identical
        Lreal = randn(LowerTriangularMatrix{NF}, spectrum)
        Lcomplex = complex.(Lreal)

        field1 = zeros(Float64, grid)
        field2 = zeros(Float64, grid)

        S = SpectralTransform(spectrum, grid; NF)

        transform!(field1, Lreal, S)
        transform!(field2, Lcomplex, S)

        for ij in eachindex(field1, field2)
            @test field1[ij] ≈ field2[ij]
        end
    end
end

@testset "Transform: NF flexibility spectral inputs" begin

    # test Float64 -> Float32
    spectrum = Spectrum(31)
    L_f32 = randn(ComplexF32, spectrum)
    L_f64 = ComplexF64.(L_f32)

    grid = OctahedralGaussianGrid(24)
    field1 = zeros(Float32, grid)
    field2 = zeros(Float32, grid)

    S = SpectralTransform(spectrum, grid, NF = Float32)

    transform!(field1, L_f64, S)    # Float64 -> Float32
    transform!(field2, L_f32, S)    # Float32 -> Float32

    for ij in eachindex(field1, field2)
        @test field1[ij] ≈ field2[ij] atol = sqrt(eps(Float32))
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

    # field3_f32 = randn(grid)
    # field3_f64 = Float64.(field3_f32)

    # transform!(L1, field3_f32, S)
    # transform!(L2, field3_f64, S)    # throws an error
end

@testset "Transform: NF flexibility for spectral outputs" begin

    # test Float32 -> Float64
    spectrum = Spectrum(31)
    grid = OctahedralGaussianGrid(24)
    field = randn(Float32, grid)

    L1 = zeros(Complex{Float32}, spectrum)
    L2 = zeros(Complex{Float64}, spectrum)

    S = SpectralTransform(spectrum, grid)

    transform!(L1, field, S)  # Float32 -> Float32
    transform!(L2, field, S)  # Float32 -> Float64

    for lm in eachindex(L1, L2)
        @test L1[lm] ≈ L2[lm]
    end

    # TODO similar to the above, NF flexibility for field outputs currently not supported because
    # the fourier transform is pre-planned and inplace with LinearAlgebra.mul! which cannot write into
    # another eltype than the FFT plan
end

@testset "Transform: Individual Legendre polynomials (inexact transforms)" begin
    @testset for trunc in spectral_resolutions_inexact
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    OctaminimalGaussianGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid,
                )

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                S = SpectralTransform(spectrum, grid; NF)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(Complex{NF}, spectrum)
                        alms[l, m] = 1

                        map = transform(alms, S)
                        alms2 = transform(map, S)

                        tol = 1.0e-3

                        for lm in eachharmonic(alms, alms2)
                            @test alms[lm] ≈ alms2[lm] atol = tol rtol = tol
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
            @testset for Grid in (
                    FullGaussianGrid,
                    FullClenshawGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                )

                # clenshaw-curtis grids are only exact for cubic truncation
                dealiasing = Grid in (FullGaussianGrid, OctahedralGaussianGrid) ? 2 : 3

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc, dealiasing))
                S = SpectralTransform(spectrum, grid; NF)

                # field that could be orography
                orography_rough = 5000 * rand(NF, grid) .^ 2
                oro_spec = transform(orography_rough, S)

                # spectral smoothing
                t = round(Int, trunc * 0.95)
                oro_spec = SpeedyTransforms.spectral_smoothing(oro_spec, 0.1, power = 1, truncation = t)

                oro_grid1 = transform(oro_spec, S)
                oro_spec1 = transform(oro_grid1, S)
                oro_grid2 = transform(oro_spec1, S)
                oro_spec2 = transform(oro_grid2, S)

                tol = NF == Float32 ? sqrt(eps(NF)) : 5.0e-7

                for lm in eachharmonic(oro_spec1, oro_spec2)
                    @test oro_spec1[lm] ≈ oro_spec2[lm] atol = tol rtol = tol
                end
                for ij in eachindex(oro_grid1, oro_grid2)
                    @test oro_grid1[ij] ≈ oro_grid2[ij] atol = tol rtol = tol
                end
            end
        end
    end
end

@testset "Transform roundtrip accuracy" begin
    # only Gaussian transforms are exact, that's why we need have different tolerenaces
    # the tolerances below were just emperically found and are kept to ensure
    # we don't introduce bugs

    tolerances = [
        1.0e-13,    # FullGaussianGrid
        1.0e-12,    # FullClenshawGrid
        1.0e-13,    # OctahedralGaussianGrid
        1.0e-12,    # OctahedralClenshawGrid
        5.0e-4,     # OctaminimalGaussianGrid
        1.0e-2,     # HEALPixGrid
        1.0e-2,     # OctaHEALPixGrid
    ]
    @testset for trunc in [42, 61]
        @testset for (i_grid, Grid) in enumerate(
                (
                    FullGaussianGrid,
                    FullClenshawGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                )
            )

            # TODO: other dealising for other grids?
            dealiasing = 3
            NF = Float64
            nlayers = 8

            spectrum = Spectrum(trunc)
            grid = Grid(SpeedyTransforms.get_nlat_half(trunc, dealiasing))
            S = SpectralTransform(spectrum, grid; NF, nlayers)

            # start in spectral space but compare in grid space to
            # avoid inaccuracies due to filtering out higher frequencies
            spec = randn(Complex{NF}, spectrum, nlayers)
            grid = transform(spec, S)

            spec_roundtrip = transform(grid, S)
            grid_roundtrip = transform(spec_roundtrip, S)

            @test grid_roundtrip ≈ grid rtol = tolerances[i_grid]
        end
    end
end
