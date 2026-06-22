@testset "SigmaCoordinates" begin

    # default and FriersonSigmaCoordinates constructors
    spectral_grid = SpectralGrid(nlayers = 8)
    @test SigmaCoordinates() isa SigmaCoordinates
    @test SigmaCoordinates(spectral_grid) isa SigmaCoordinates
    @test FriersonSigmaCoordinates(spectral_grid) isa SigmaCoordinates

    # constructor from vector and from range
    @test SigmaCoordinates([0, 0.4, 0.6, 1.0]) isa SigmaCoordinates
    @test SigmaCoordinates(0:0.25:1) isa SigmaCoordinates

    # get_nlayers
    σ = SigmaCoordinates(spectral_grid)
    @test SpeedyWeather.get_nlayers(σ) == 8

    # get_σ_half returns the stored half levels
    σ3 = SigmaCoordinates([0.0, 0.4, 0.6, 1.0])
    @test SpeedyWeather.get_σ_half(σ3) ≈ Float32[0, 0.4, 0.6, 1]

    # boundary conditions: top = 0, surface = 1, strictly increasing
    for nlayers in (3, 5, 8)
        σ = SigmaCoordinates(SpectralGrid(nlayers = nlayers))
        σh = Array(σ.σ_half)
        @test σh[1] == 0
        @test σh[end] == 1
        @test all(diff(σh) .> 0)
    end

    # σ_full = midpoints of σ_half
    spectral_grid4 = SpectralGrid(nlayers = 4)
    σ4 = SigmaCoordinates(spectral_grid4)
    σh = Array(σ4.σ_half)
    @test Array(σ4.σ_full) ≈ (σh[1:end-1] + σh[2:end]) / 2

    # σ_thickness sums to 1 (covers full atmosphere)
    @test sum(σ4.σ_thickness) ≈ 1

    # pressure and pressure_thickness
    p_surf = 1.0e5
    for k in 1:4
        @test SpeedyWeather.pressure(k, p_surf, σ4) ≈ σ4.σ_full[k] * p_surf
        @test SpeedyWeather.pressure_thickness(k, p_surf, σ4) ≈ σ4.σ_thickness[k] * p_surf
    end

    # sigma_half_spacing: default is linear
    σh_linear = SpeedyWeather.sigma_half_spacing(4)
    @test σh_linear ≈ collect(0:0.25:1)

    # sigma_half_spacing with custom profile
    σh_custom = SpeedyWeather.sigma_half_spacing(4, x -> x^2)
    @test σh_custom[1] == 0
    @test σh_custom[end] ≈ 1

    # Frierson levels satisfy boundary conditions
    σ_frier = FriersonSigmaCoordinates(spectral_grid4)
    @test Array(σ_frier.σ_half)[1] == 0
    @test Array(σ_frier.σ_half)[end] ≈ 1
    @test SpeedyWeather.get_nlayers(σ_frier) == 4

    # sigma_okay rejects invalid inputs
    @test_throws AssertionError SpeedyWeather.sigma_okay(3, [-0.1, 0.4, 0.6, 1.0]) # first < 0
    @test_throws AssertionError SpeedyWeather.sigma_okay(3, [0.0, 0.4, 0.6, 0.9])  # last != 1
    @test_throws AssertionError SpeedyWeather.sigma_okay(2, [0.0, 0.4, 0.6, 1.0])  # wrong nlayers
    @test_throws AssertionError SpeedyWeather.sigma_okay(3, [0.0, 0.6, 0.4, 1.0])  # not increasing

    # automatic levels in Geometry
    G = Geometry(spectral_grid4)
    @test length(G.σ_levels_half) == 5
    @test length(G.σ_levels_full) == 4

    # manual levels from vector, infer nlayers
    σ_manual = SigmaCoordinates([0, 0.4, 0.6, 1.0])
    spectral_grid3 = SpectralGrid(nlayers = SpeedyWeather.get_nlayers(σ_manual))
    G = Geometry(spectral_grid3, vertical_coordinates = σ_manual)
    @test spectral_grid3.nlayers == 3
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3

    # specify SpectralGrid and levels together
    σ_both = SigmaCoordinates(spectral_grid3, [0, 0.4, 0.6, 1.0])
    G = Geometry(spectral_grid3, vertical_coordinates = σ_both)
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3
end

@testset "pressure and pressure_thickness interface" begin

    nlayers = 4
    σ_half_vals = Float32[0, 0.2, 0.5, 0.8, 1.0]
    p_surf = 1.0f5
    spectral_grid = SpectralGrid(nlayers = nlayers)

    σ = SigmaCoordinates(spectral_grid, σ_half_vals)

    # SigmaPressureCoordinates with transition = 1 everywhere => A = 0, B = σ (pure sigma)
    S = SigmaPressureCoordinates(spectral_grid; σ_half = σ_half_vals, transition = _ -> 1.0)
    @test all(Array(S.A_half) .≈ 0)
    @test all(Array(S.A_full) .≈ 0)

    # pressure and pressure_thickness must agree between both coordinate types
    for k in 1:nlayers
        @test SpeedyWeather.pressure(k, p_surf, σ) ≈ SpeedyWeather.pressure(k, p_surf, S)
        @test SpeedyWeather.pressure_thickness(k, p_surf, σ) ≈ SpeedyWeather.pressure_thickness(k, p_surf, S)
    end
end

@testset "SigmaPressureCoordinates" begin

    # default construction
    spectral_grid = SpectralGrid(nlayers = 4)
    S = SigmaPressureCoordinates(spectral_grid)
    @test S isa SigmaPressureCoordinates
    @test SpeedyWeather.get_nlayers(S) == 4

    # A + B = σ at half and full levels (split is exact)
    σ_half = SpeedyWeather.sigma_half_spacing(spectral_grid.nlayers)
    @test Array(S.A_half) + Array(S.B_half) ≈ σ_half
    σ_full = (σ_half[2:end] + σ_half[1:(end-1)]) / 2
    @test Array(S.A_full) + Array(S.B_full) ≈ σ_full

    # custom reference pressure
    S_p = SigmaPressureCoordinates(spectral_grid; reference_pressure = 5.0e4)
    @test S_p.reference_pressure ≈ 5.0f4

    # custom σ_half levels
    spectral_grid3 = SpectralGrid(nlayers = 3)
    S3 = SigmaPressureCoordinates(spectral_grid3; σ_half = [0.0, 0.3, 0.7, 1.0])
    @test SpeedyWeather.get_nlayers(S3) == 3
    @test length(S3.A_half) == 4
    @test length(S3.B_full) == 3

    # pure sigma: transition = σ -> 1  =>  A = 0, B = σ
    S_sigma = SigmaPressureCoordinates(spectral_grid; transition = σ -> 1.0)
    @test all(Array(S_sigma.A_half) .≈ 0)
    @test all(Array(S_sigma.A_full) .≈ 0)
    @test Array(S_sigma.B_half) ≈ σ_half
    @test Array(S_sigma.B_full) ≈ σ_full

    # pure pressure: transition = σ -> 0  =>  A = σ, B = 0
    S_pres = SigmaPressureCoordinates(spectral_grid; transition = σ -> 0.0)
    @test all(Array(S_pres.B_half) .≈ 0)
    @test all(Array(S_pres.B_full) .≈ 0)
    @test Array(S_pres.A_half) ≈ σ_half
    @test Array(S_pres.A_full) ≈ σ_full

    # use inside Geometry
    G = Geometry(spectral_grid, vertical_coordinates = S)
    @test G isa Geometry
    @test length(G.σ_levels_half) == 5
    @test length(G.σ_levels_full) == 4

    # Geometry constructed directly with SigmaPressureCoordinates and custom levels
    G3 = Geometry(spectral_grid3, vertical_coordinates = S3)
    @test G3 isa Geometry
    @test length(G3.σ_levels_half) == 4
    @test length(G3.σ_levels_full) == 3
end
