@testset "SigmaCoordinates" begin

    spectral_grid = SpectralGrid(nlayers = 8)

    # no-SpectralGrid constructors
    @test SigmaCoordinates() isa SigmaCoordinates
    @test SigmaCoordinates([0, 0.4, 0.6, 1.0]) isa SigmaCoordinates
    @test SigmaCoordinates(0:0.25:1) isa SigmaCoordinates

    # SpectralGrid first, vector second
    @test SigmaCoordinates(spectral_grid) isa SigmaCoordinates
    @test SigmaCoordinates(spectral_grid, collect(0:0.125:1)) isa SigmaCoordinates

    # SpectralGrid first, function second
    @test SigmaCoordinates(spectral_grid, x -> x^2) isa SigmaCoordinates
    @test SigmaCoordinates(spectral_grid, SpeedyWeather.frierson_profile) isa SigmaCoordinates

    # get_nlayers
    σ = SigmaCoordinates(spectral_grid)
    @test SpeedyWeather.get_nlayers(σ) == 8

    # get_σ_half, get_σ_full, get_σ_thickness
    σ3 = SigmaCoordinates([0.0, 0.4, 0.6, 1.0])
    @test SpeedyWeather.get_σ_half(σ3) ≈ Float32[0, 0.4, 0.6, 1]
    @test SpeedyWeather.get_σ_full(σ3) ≈ Float32[0.2, 0.5, 0.8]
    @test SpeedyWeather.get_σ_thickness(σ3) ≈ Float32[0.4, 0.2, 0.4]
    @test sum(SpeedyWeather.get_σ_thickness(σ3)) ≈ 1

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

@testset "FriersonSigmaCoordinates" begin

    spectral_grid = SpectralGrid(nlayers = 4)

    # no-SpectralGrid and SpectralGrid constructors
    @test FriersonSigmaCoordinates() isa SigmaCoordinates
    @test FriersonSigmaCoordinates(spectral_grid) isa SigmaCoordinates

    # boundary conditions
    σ_frier = FriersonSigmaCoordinates(spectral_grid)
    @test Array(σ_frier.σ_half)[1] == 0
    @test Array(σ_frier.σ_half)[end] ≈ 1
    @test SpeedyWeather.get_nlayers(σ_frier) == 4
end

@testset "SigmaPressureCoordinates" begin

    spectral_grid = SpectralGrid(nlayers = 4)

    # no-SpectralGrid constructor
    @test SigmaPressureCoordinates() isa SigmaPressureCoordinates

    # SpectralGrid only (default σ_half, default transition)
    S = SigmaPressureCoordinates(spectral_grid)
    @test S isa SigmaPressureCoordinates
    @test SpeedyWeather.get_nlayers(S) == 4

    # get_σ_full and get_σ_thickness return the sigma-terrain-following (B) components
    @test SpeedyWeather.get_σ_full(S) ≈ Array(S.B_full)
    @test SpeedyWeather.get_σ_thickness(S) ≈ Array(S.B_thickness)

    # with pure sigma (transition=1), get_σ_full and get_σ_thickness match SigmaCoordinates
    S_sigma = SigmaPressureCoordinates(spectral_grid; transition = _ -> 1.0)
    σ_sigma = SigmaCoordinates(spectral_grid)
    @test SpeedyWeather.get_σ_full(S_sigma) ≈ Array(SpeedyWeather.get_σ_full(σ_sigma))
    @test SpeedyWeather.get_σ_thickness(S_sigma) ≈ Array(SpeedyWeather.get_σ_thickness(σ_sigma))

    # SpectralGrid + vector (second positional)
    S3 = SigmaPressureCoordinates(SpectralGrid(nlayers = 3), [0.0, 0.3, 0.7, 1.0])
    @test SpeedyWeather.get_nlayers(S3) == 3
    @test length(S3.A_half) == 4
    @test length(S3.B_full) == 3

    # SpectralGrid + function (second positional)
    S_fn = SigmaPressureCoordinates(spectral_grid, SpeedyWeather.frierson_profile)
    @test S_fn isa SigmaPressureCoordinates
    @test SpeedyWeather.get_nlayers(S_fn) == 4

    # keyword args: reference_pressure
    S_p = SigmaPressureCoordinates(spectral_grid; reference_pressure = 5.0e4)
    @test S_p.reference_pressure ≈ 5.0f4

    # keyword args: transition
    # pure sigma: transition = σ -> 1  =>  A = 0, B = σ
    σ_half = SpeedyWeather.sigma_half_spacing(spectral_grid.nlayers)
    σ_full = (σ_half[2:end] + σ_half[1:(end-1)]) / 2
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

    # A + B = σ at half and full levels (split is exact)
    @test Array(S.A_half) + Array(S.B_half) ≈ σ_half
    @test Array(S.A_full) + Array(S.B_full) ≈ σ_full

    # use inside Geometry
    G = Geometry(spectral_grid, vertical_coordinates = S)
    @test G isa Geometry
    @test length(G.σ_levels_half) == 5
    @test length(G.σ_levels_full) == 4

    G3 = Geometry(SpectralGrid(nlayers = 3), vertical_coordinates = S3)
    @test G3 isa Geometry
    @test length(G3.σ_levels_half) == 4
    @test length(G3.σ_levels_full) == 3
end

@testset "CubicSigmaPressureCoordinates" begin

    spectral_grid = SpectralGrid(nlayers = 8)
    pressure_only_above = 0.2f0
    σ_only_below = 0.8f0

    # no-SpectralGrid constructor
    @test CubicSigmaPressureCoordinates() isa SigmaPressureCoordinates

    # SpectralGrid only (uses default thresholds)
    H = CubicSigmaPressureCoordinates(spectral_grid; pressure_only_above, σ_only_below)
    @test H isa SigmaPressureCoordinates
    @test SpeedyWeather.get_nlayers(H) == 8

    # SpectralGrid + vector (second positional)
    H_vec = CubicSigmaPressureCoordinates(spectral_grid, collect(range(0.0, 1.0, 9)); pressure_only_above, σ_only_below)
    @test H_vec isa SigmaPressureCoordinates
    @test SpeedyWeather.get_nlayers(H_vec) == 8

    # SpectralGrid + function (second positional)
    H_fn = CubicSigmaPressureCoordinates(spectral_grid, SpeedyWeather.frierson_profile; pressure_only_above, σ_only_below)
    @test H_fn isa SigmaPressureCoordinates

    # cubic_transition: pure pressure at top, pure sigma at bottom
    @test SpeedyWeather.cubic_transition(0.0f0; pressure_only_above, σ_only_below) == 0.0f0
    @test SpeedyWeather.cubic_transition(pressure_only_above; pressure_only_above, σ_only_below) == 0.0f0
    @test SpeedyWeather.cubic_transition(σ_only_below; pressure_only_above, σ_only_below) == 1.0f0
    @test SpeedyWeather.cubic_transition(1.0f0; pressure_only_above, σ_only_below) == 1.0f0

    # cubic_transition: smooth and monotone in between
    σ_mid = range(pressure_only_above, σ_only_below, 20)
    t_mid = SpeedyWeather.cubic_transition.(σ_mid; pressure_only_above, σ_only_below)
    @test all(0 .<= t_mid .<= 1)
    @test all(diff(t_mid) .>= 0)

    # C¹ continuity: near-zero derivative at both thresholds
    ε = 1.0f-4
    @test SpeedyWeather.cubic_transition(pressure_only_above + ε; pressure_only_above, σ_only_below) ≈ 0 atol = 1e-6
    @test SpeedyWeather.cubic_transition(σ_only_below - ε; pressure_only_above, σ_only_below) ≈ 1 atol = 1e-2

    # thresholds respected in the coordinate itself
    B_half = Array(H.B_half)
    A_half = Array(H.A_half)
    σ_half = A_half + B_half
    @test all(B_half[σ_half .<= pressure_only_above] .≈ 0)   # pure pressure above pressure_only_above
    @test all(A_half[σ_half .>= σ_only_below] .≈ 0)          # pure sigma below σ_only_below

    # custom thresholds keyword args
    H2 = CubicSigmaPressureCoordinates(spectral_grid; pressure_only_above = 0.1, σ_only_below = 0.9)
    B2 = Array(H2.B_half)
    A2 = Array(H2.A_half)
    σ2 = A2 + B2
    @test all(B2[σ2 .<= 0.1] .≈ 0)
    @test all(A2[σ2 .>= 0.9] .≈ 0)

    # smoothstep coefficients are independent of thresholds: midpoint of transition always = 0.5
    @test SpeedyWeather.cubic_transition(0.5f0; pressure_only_above = 0.0f0, σ_only_below = 1.0f0) ≈ 0.5f0
    @test SpeedyWeather.cubic_transition(0.5f0; pressure_only_above = 0.3f0, σ_only_below = 0.7f0) ≈ 0.5f0

    # works inside Geometry
    G = Geometry(spectral_grid, vertical_coordinates = H)
    @test G isa Geometry
    @test length(G.σ_levels_half) == 9
end

@testset "pressure and pressure_thickness interface" begin

    nlayers = 4
    σ_half_vals = Float32[0, 0.2, 0.5, 0.8, 1.0]
    p_surf = 1.0f5
    spectral_grid = SpectralGrid(nlayers = nlayers)

    σ = SigmaCoordinates(spectral_grid, σ_half_vals)

    # SigmaPressureCoordinates with transition = 1 everywhere => A = 0, B = σ (pure sigma)
    S = SigmaPressureCoordinates(spectral_grid, σ_half_vals; transition = _ -> 1.0)
    @test all(Array(S.A_half) .≈ 0)
    @test all(Array(S.A_full) .≈ 0)

    # pressure and pressure_thickness must agree between both coordinate types
    for k in 1:nlayers
        @test SpeedyWeather.pressure(k, p_surf, σ) ≈ SpeedyWeather.pressure(k, p_surf, S)
        @test SpeedyWeather.pressure_thickness(k, p_surf, σ) ≈ SpeedyWeather.pressure_thickness(k, p_surf, S)
    end
end

@testset "sigma function" begin

    nlayers = 4
    σ_half_vals = Float32[0, 0.2, 0.5, 0.8, 1.0]
    spectral_grid = SpectralGrid(nlayers = nlayers)

    σ = SigmaCoordinates(spectral_grid, σ_half_vals)

    # sigma returns σ_full for SigmaCoordinates
    for k in 1:nlayers
        @test SpeedyWeather.sigma(k, σ) ≈ σ.σ_full[k]
    end

    # sigma is surface-pressure independent (unlike pressure)
    S_default = SigmaPressureCoordinates(spectral_grid, σ_half_vals)
    for k in 1:nlayers
        @test SpeedyWeather.sigma(k, S_default) ≈ S_default.A_full[k] + S_default.B_full[k]
    end

    # sigma agrees between SigmaCoordinates and SigmaPressureCoordinates for all transitions:
    # A + B = σ always, so sigma is the same nominal level regardless of how it is split
    for transition in (_ -> 0.0, _ -> 0.5, _ -> 1.0, σ -> σ)
        S = SigmaPressureCoordinates(spectral_grid, σ_half_vals; transition)
        for k in 1:nlayers
            @test SpeedyWeather.sigma(k, σ) ≈ SpeedyWeather.sigma(k, S)
        end
    end

    # sigma lies in (0, 1] for all interior and surface layers
    for k in 1:nlayers
        @test 0 < SpeedyWeather.sigma(k, σ) <= 1
    end
end

@testset "pressure_half, pressure_above, pressure_below" begin

    nlayers = 4
    spectral_grid = SpectralGrid(nlayers = nlayers)
    p_surf = 1.0f5

    for coord in (SigmaCoordinates(spectral_grid), SigmaPressureCoordinates(spectral_grid))

        # pressure_half: nlayers+1 interfaces, increasing from top (0 Pa) to surface
        p_half = [SpeedyWeather.pressure_half(k, p_surf, coord) for k in 1:(nlayers + 1)]
        @test p_half[1] == 0                   # top of atmosphere
        @test p_half[end] ≈ p_surf             # surface
        @test all(diff(p_half) .> 0)           # strictly increasing downward

        for k in 1:nlayers
            # pressure_above / pressure_below bracket the full level
            p_above = SpeedyWeather.pressure_above(k, p_surf, coord)
            p_full  = SpeedyWeather.pressure(k, p_surf, coord)
            p_below = SpeedyWeather.pressure_below(k, p_surf, coord)
            @test p_above <= p_full <= p_below

            # pressure_thickness equals the difference between the two interfaces
            Δp = SpeedyWeather.pressure_thickness(k, p_surf, coord)
            @test Δp ≈ p_below - p_above

            # pressure_above and pressure_below are consistent with pressure_half
            @test p_above ≈ SpeedyWeather.pressure_half(k,     p_surf, coord)
            @test p_below ≈ SpeedyWeather.pressure_half(k + 1, p_surf, coord)
        end

        # thicknesses sum to surface pressure
        Δp_sum = sum(SpeedyWeather.pressure_thickness(k, p_surf, coord) for k in 1:nlayers)
        @test Δp_sum ≈ p_surf
    end
end

@testset "Reference pressure mismatch warning" begin
    spectral_grid = SpectralGrid(nlayers = 4)

    # mismatched: coordinate has 5e4 Pa, atmosphere defaults to 1e5 Pa
    S = SigmaPressureCoordinates(spectral_grid; reference_pressure = 5.0e4)
    model = PrimitiveDryModel(spectral_grid; geometry = Geometry(spectral_grid; vertical_coordinates = S))
    @test_logs (:warn, r"Reference pressure") initialize!(model)

    # matched: no warning expected
    S_match = SigmaPressureCoordinates(spectral_grid; reference_pressure = 1.0e5)
    model_match = PrimitiveDryModel(spectral_grid; geometry = Geometry(spectral_grid; vertical_coordinates = S_match))
    @test_logs initialize!(model_match)
end

@testset "Held-Suarez with varying layers and sigma spacing" begin
    for nlayers in (16, 24)     # don't test more layers regularly on CI but 64 should work too
        spectral_grid = SpectralGrid(trunc = 31; nlayers)

        σ_half_configs = (
            ("default",   SpeedyWeather.sigma_half_spacing(nlayers)),
            ("Frierson",  SpeedyWeather.sigma_half_spacing(nlayers, SpeedyWeather.frierson_profile)),
        )

        for (label, σ_half) in σ_half_configs
            coords = (
                ("SigmaCoordinates/$label",          SigmaCoordinates(spectral_grid, σ_half)),
                # like sigma only for now with A=0
                ("SigmaPressureCoordinates/$label",  SigmaPressureCoordinates(spectral_grid, σ_half; transition = _ -> 1.0)),
            )

            for (name, coord) in coords
                @testset "$name (nlayers=$nlayers)" begin
                    @info "$name (nlayers=$nlayers)"
                    model = PrimitiveDryModel(
                        spectral_grid;
                        geometry  = Geometry(spectral_grid; vertical_coordinates = coord),
                        forcing   = HeldSuarez(spectral_grid),
                        drag      = LinearDrag(spectral_grid),
                        dynamics_only = true,
                        orography = EarthOrography(spectral_grid),
                    )
                    simulation = initialize!(model)
                    run!(simulation; period = Day(5))               # don't test longer on CI but otherwise could be Day(20) too
                    if label == "default" && nlayers > 32           # current this blows up but that might be numerically expected
                        @test_broken model.feedback.nans_detected == false
                    else
                        @test model.feedback.nans_detected == false
                    end
                end
            end
        end
    end
end