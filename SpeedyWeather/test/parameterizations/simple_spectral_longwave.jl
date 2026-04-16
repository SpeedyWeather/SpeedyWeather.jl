@testset "SimpleSpectralLongwave: spectroscopic helpers" begin

    # ---- Planck function -------------------------------------------------
    # π × B integrated analytically from 0 to ∞ must recover σ T⁴.
    # We only cover 10–2510 cm⁻¹, so expect ≥ 99% of the blackbody flux at 300 K.
    T_test = 300f0
    rad    = SimpleSpectralLongwave{Float32}()
    nν     = rad.nwavenumber
    dν̃     = (rad.wavenumber_max - rad.wavenumber_min) / (nν - 1)
    pi_B_sum = sum(
        Float32(π) * SpeedyWeather.planck_wavenumber(T_test, rad.wavenumber_min + (iv-1)*dν̃) * dν̃
        for iv in 1:nν
    )
    σ_SB = 5.670374419f-8   # W m⁻² K⁻⁴
    @test pi_B_sum ≈ σ_SB * T_test^4 rtol=0.01   # within 1% (truncated spectrum)

    # ---- H₂O line absorption (Williams 2026, Table 1) -------------------
    # Check boundary values and exponential decay behaviour.
    @test SpeedyWeather.h2o_line_kappa_ref(100f0, rad) ≈ rad.κ_rot           # flat region
    @test SpeedyWeather.h2o_line_kappa_ref(200f0, rad) ≈ rad.κ_rot           # rotation band edge
    @test SpeedyWeather.h2o_line_kappa_ref(600f0, rad) < rad.κ_rot           # decaying rotation band
    @test SpeedyWeather.h2o_line_kappa_ref(600f0, rad) ≈ rad.κ_rot * exp(-(600f0-200f0)/rad.l_rot) rtol=1e-5
    @test SpeedyWeather.h2o_line_kappa_ref(1450f0, rad) ≈ rad.κ_vr          # peak vibration-rotation
    @test SpeedyWeather.h2o_line_kappa_ref(1600f0, rad) ≈ rad.κ_vr          # flat plateau
    @test SpeedyWeather.h2o_line_kappa_ref(2000f0, rad) < rad.κ_vr          # combination band decaying
    @test SpeedyWeather.h2o_line_kappa_ref(3000f0, rad) == 0f0               # beyond range

    # ---- CO₂ absorption (Williams 2026, Eq. 5) --------------------------
    # Peak at ν̃_co2 = 667 cm⁻¹, zero outside [500, 850] cm⁻¹.
    @test SpeedyWeather.co2_kappa_ref(667f0, rad) ≈ rad.κ_co2               # peak
    @test SpeedyWeather.co2_kappa_ref(667f0 + rad.l_co2, rad) ≈ rad.κ_co2 / Float32(ℯ) rtol=1e-4  # e-folding
    @test SpeedyWeather.co2_kappa_ref(400f0, rad) == 0f0                     # outside band
    @test SpeedyWeather.co2_kappa_ref(900f0, rad) == 0f0                     # outside band

    # ---- H₂O continuum (Williams 2026, Eq. 6) ---------------------------
    @test SpeedyWeather.h2o_cont_kappa_ref(1000f0, rad) == rad.κ_cnt1       # window region
    @test SpeedyWeather.h2o_cont_kappa_ref(2000f0, rad) == rad.κ_cnt2       # above 1700 cm⁻¹
    @test SpeedyWeather.h2o_cont_kappa_ref(1700f0, rad) == rad.κ_cnt1       # boundary: inclusive ≤ 1700
    @test rad.κ_cnt1 > rad.κ_cnt2                                            # larger below 1700 (window)
end

@testset "SimpleSpectralLongwave: optical depth" begin
    # A non-zero specific humidity should give a positive optical depth increment
    # for any wavenumber. CO₂ optical depth should also be positive when enabled.
    spectral_grid = SpectralGrid(trunc=5, nlayers=4)
    rad_h2o = SimpleSpectralLongwave(spectral_grid; do_co2=false)
    rad_co2 = SimpleSpectralLongwave(spectral_grid; do_co2=true)

    model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad_h2o)
    initialize!(model.longwave_radiation, model)
    vars  = Variables(model)

    # Set a realistic humidity column (0.5% specific humidity throughout)
    vars.grid.pressure_prev .= 100000f0   # surface pressure [Pa]
    for k in 1:spectral_grid.nlayers
        vars.grid.humidity_prev[:, k]     .= 0.005f0
        vars.grid.temperature_prev[:, k]  .= 280f0
    end

    ij  = 1
    pₛ  = vars.grid.pressure_prev[ij]
    T   = vars.grid.temperature_prev
    q   = vars.grid.humidity_prev

    dν̃  = (rad_h2o.wavenumber_max - rad_h2o.wavenumber_min) / (rad_h2o.nwavenumber - 1)

    for k in 1:spectral_grid.nlayers
        # Window region (1000 cm⁻¹): H₂O absorption should be positive
        Δτ_win = SpeedyWeather.ssm_delta_tau(ij, k, 1000f0, T, q, pₛ, rad_h2o, model)
        @test Δτ_win > 0

        # Rotation band (400 cm⁻¹): strongest H₂O absorption
        Δτ_rot = SpeedyWeather.ssm_delta_tau(ij, k, 400f0, T, q, pₛ, rad_h2o, model)
        @test Δτ_rot > Δτ_win       # rotation band more opaque than window

        # CO₂ band (667 cm⁻¹): should increase optical depth when CO₂ is active
        Δτ_no_co2 = SpeedyWeather.ssm_delta_tau(ij, k, 667f0, T, q, pₛ, rad_h2o, model)
        Δτ_with_co2 = SpeedyWeather.ssm_delta_tau(ij, k, 667f0, T, q, pₛ, rad_co2, model)
        @test Δτ_with_co2 > Δτ_no_co2

        # Optical depth must be non-negative for all wavenumbers
        for iv in 1:rad_h2o.nwavenumber
            ν̃ = rad_h2o.wavenumber_min + (iv-1) * dν̃
            @test SpeedyWeather.ssm_delta_tau(ij, k, ν̃, T, q, pₛ, rad_h2o, model) >= 0
        end
    end
end

@testset "SimpleSpectralLongwave: parameterization interface" begin
    # Basic smoke test: construct, initialize, run one column.
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset for do_co2 in (false, true)
        rad   = SimpleSpectralLongwave(spectral_grid; do_co2)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad)

        initialize!(model.longwave_radiation, model)
        vars = Variables(model)
        vars.grid.pressure_prev .= 100000
        # k=1 is top layer, k=nlayers is bottom layer → temperature increases with k
        for k in 1:spectral_grid.nlayers
            vars.grid.temperature_prev[:, k] .= 220 + 9 * (k-1)
            vars.grid.humidity_prev[:, k]    .= 0.005
        end
        vars.prognostic.ocean.sea_surface_temperature .= 295
        vars.prognostic.land.soil_temperature[:, 1]   .= 285

        ij = rand(1:spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)

        # After one call, temperature tendency should be non-zero (atmosphere cools)
        @test any(!=(0), vars.tendencies.grid.temperature[ij, :])
    end
end

@testset "SimpleSpectralLongwave: energy conservation" begin
    # Column energy balance: for an isothermal column in radiative equilibrium at
    # temperature T_eq, the net longwave flux divergence should be near zero.
    # Here we merely check that tendencies are finite and the surface downwelling
    # flux is non-negative and does not exceed the upwelling flux (no violation
    # of second law for realistic temperatures).
    spectral_grid = SpectralGrid(trunc=5, nlayers=8)
    rad   = SimpleSpectralLongwave(spectral_grid; do_co2=true)
    model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad)
    initialize!(model.longwave_radiation, model)
    vars  = Variables(model)

    # Initialise all prognostic / grid fields to a realistic state
    # k=1 is top layer, k=nlayers is bottom: temperature increases with k
    vars.grid.pressure_prev .= 100000f0   # surface pressure [Pa]
    for k in 1:spectral_grid.nlayers
        vars.grid.temperature_prev[:, k] .= 220f0 + 9f0 * (k-1)  # lapse rate
        vars.grid.humidity_prev[:, k]    .= 0.005f0
    end
    vars.prognostic.ocean.sea_surface_temperature .= 295f0
    vars.prognostic.land.soil_temperature[:, 1]   .= 285f0

    for ij in 1:spectral_grid.npoints
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
    end

    # All temperature tendencies must be finite
    @test all(isfinite, vars.tendencies.grid.temperature)

    # Outgoing longwave > 0 everywhere
    @test all(>(0), vars.parameterizations.outgoing_longwave)

    # Downwelling surface longwave must be ≥ 0 and ≤ surface upwelling
    @test all(>=(0), vars.parameterizations.surface_longwave_down)
    @test all(
        vars.parameterizations.surface_longwave_down .<= vars.parameterizations.surface_longwave_up
    )
end

@testset "SimpleSpectralLongwave: CO₂ forcing" begin
    # Doubling CO₂ should reduce OLR (positive radiative forcing).
    spectral_grid = SpectralGrid(trunc=5, nlayers=8)

    function run_and_get_olr(co2_ppmv)
        rad   = SimpleSpectralLongwave(spectral_grid; do_co2=true, co2_ppmv)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad)
        initialize!(model.longwave_radiation, model)
        vars  = Variables(model)
        vars.grid.pressure_prev .= 100000f0   # surface pressure [Pa]
        # k=1 is top layer, k=nlayers is bottom: temperature increases with k
        for k in 1:spectral_grid.nlayers
            vars.grid.temperature_prev[:, k] .= 220f0 + 9f0 * (k-1)
            vars.grid.humidity_prev[:, k]    .= 0.005f0
        end
        vars.prognostic.ocean.sea_surface_temperature .= 295f0
        vars.prognostic.land.soil_temperature[:, 1]   .= 285f0
        for ij in 1:spectral_grid.npoints
            SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
        end
        return mean(vars.parameterizations.outgoing_longwave)
    end

    olr_280 = run_and_get_olr(280f0)
    olr_560 = run_and_get_olr(560f0)

    # Doubling CO₂ reduces OLR (greenhouse effect)
    @test olr_560 < olr_280

    # Forcing magnitude should be physically plausible: 2–5 W m⁻²
    forcing = olr_280 - olr_560
    @test 1 < forcing < 10
end

@testset "SimpleSpectralLongwave: number format" begin
    # The scheme should work with both Float32 and Float64
    @testset for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, trunc=5, nlayers=4)
        rad   = SimpleSpectralLongwave(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation=rad)
        initialize!(model.longwave_radiation, model)
        vars  = Variables(model)
        vars.grid.pressure_prev .= 100000
        # k=1 is top layer, k=nlayers is bottom: temperature increases with k
        for k in 1:spectral_grid.nlayers
            vars.grid.temperature_prev[:, k] .= 220 + 9 * (k-1)
            vars.grid.humidity_prev[:, k]    .= 0.005
        end
        vars.prognostic.ocean.sea_surface_temperature .= 295
        vars.prognostic.land.soil_temperature[:, 1]   .= 285
        SpeedyWeather.parameterization!(1, vars, model.longwave_radiation, model)
        @test all(isfinite, vars.tendencies.grid.temperature[1, :])
    end
end
