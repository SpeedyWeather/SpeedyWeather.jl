# Testing and verification items 4-6 of
# docs/dev/2026-08/hybrid-sigma-pressure-coordinates-part-2.md: the generalisation of the
# dynamical core (vertical_integration!, vertical_velocity!, vertical_advection!,
# temperature_grid_tendency!, surface_pressure_grid_tendency!) to hybrid sigma-pressure
# coordinates must reduce exactly (up to the extra log/exp the hybrid path evaluates) to the
# existing SigmaCoordinates results in the pure-sigma limit (transition = _ -> 1), mass must be
# conserved by the coordinate weights, and a genuine (non-degenerate) hybrid coordinate must
# not blow up in a real integration (the regression test for the NaN bug fixed in this PR).

@testset "hybrid vs sigma: identical dynamics_tendencies! in the pure-sigma limit" begin
    # item 4: build a SigmaCoordinates model and a SigmaPressureCoordinates(transition = _ ->
    # 1.0) model (the pure-sigma limit of the hybrid code path), run the SigmaCoordinates model
    # forward briefly to get a physically sensible (non-trivial) atmospheric state, then copy
    # that single state into both leapfrog steps of both models, and compare one call of
    # dynamics_tendencies! field by field.
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)

    function make_model(vertical_coordinates)
        geometry = isnothing(vertical_coordinates) ? Geometry(spectral_grid) : Geometry(spectral_grid; vertical_coordinates)
        model = PrimitiveDryModel(
            spectral_grid;
            geometry,
            orography = EarthOrography(spectral_grid),
            forcing = HeldSuarez(spectral_grid),
            drag = LinearDrag(spectral_grid),
            dynamics_only = true,
        )
        model.feedback.verbose = false
        return model
    end

    sigma_model = make_model(nothing)
    sigma_sim = initialize!(sigma_model)
    run!(sigma_sim; period = Hour(1))    # `run!` leaves the prognostic state unscaled (scale = 1)

    coord = SigmaPressureCoordinates(spectral_grid; transition = _ -> 1.0)
    hybrid_model = make_model(coord)
    hybrid_sim = initialize!(hybrid_model)

    # `implicit.temp_profile` is a linearisation reference computed once in `first_timesteps!`
    # from the *running* average temperature (`reinitialize!`); since `hybrid_model` is only
    # `initialize!`d here (never stepped through `run!`), it would otherwise stay at its
    # zero-allocated default while `sigma_model`'s reflects the evolved state above. This
    # reference profile is not itself part of the coordinate-dependent physics under test
    # (it is subtracted from and added back to temperature symmetrically), so syncing it
    # avoids polluting the comparison below with an unrelated, incidental mismatch.
    hybrid_model.implicit.temp_profile .= sigma_model.implicit.temp_profile
    hybrid_model.implicit.Δt[] = sigma_model.implicit.Δt[]

    # copy the (single, current-step) evolved state of the sigma model into both leapfrog
    # steps of both models, so that dynamics_tendencies! starts from an identical state that
    # does not depend on either model's own leapfrog history
    for sim in (sigma_sim, hybrid_sim)
        vars = sim.variables
        for name in (:vorticity, :divergence, :temperature, :pressure)
            src = getfield(sigma_sim.variables.prognostic, name)
            dst = getfield(vars.prognostic, name)
            SpeedyWeather.get_step(dst, 2) .= SpeedyWeather.get_step(src, 2)
            SpeedyWeather.get_step(dst, 1) .= SpeedyWeather.get_step(src, 2)
        end
    end

    for (sim, model) in ((sigma_sim, sigma_model), (hybrid_sim, hybrid_model))
        vars = sim.variables
        # the dynamical core assumes vorticity/divergence are radius-scaled throughout (as
        # they are while `run!`'s time loop is executing; `run!` only unscales once at the
        # very end for user-facing output) -- redo that scaling before calling
        # dynamics_tendencies! directly, mirroring what `run!` does before its time loop
        SpeedyWeather.scale_prognostic!(vars, model.planet.radius)
        # two consecutive transforms make grid step 1 == grid step 2 == the copied snapshot
        # (the first call moves the *old* grid state into step 1, the second moves the
        # now-consistent step 2 into step 1) -- this removes any dependence of the comparison
        # below on the two models' distinct leapfrog history
        transform!(vars, model)
        transform!(vars, model)
        SpeedyWeather.reset_tendencies!(vars, model.time_stepping)
    end

    SpeedyWeather.dynamics_tendencies!(sigma_sim.variables, sigma_model)
    SpeedyWeather.dynamics_tendencies!(hybrid_sim.variables, hybrid_model)

    vars_σ = sigma_sim.variables
    vars_h = hybrid_sim.variables
    rtol = 1.0f-4

    # vertical_integration!: u_mean_grid, v_mean_grid ΔB-weighted; div_mean_grid, div_sum_above
    # δ-weighted; div_mean (spectral) stays Δσ-weighted; div_mean_correction is the difference
    @test Array(vars_h.dynamics.u_mean_grid) ≈ Array(vars_σ.dynamics.u_mean_grid) rtol = rtol
    @test Array(vars_h.dynamics.v_mean_grid) ≈ Array(vars_σ.dynamics.v_mean_grid) rtol = rtol
    @test Array(vars_h.dynamics.div_mean_grid) ≈ Array(vars_σ.dynamics.div_mean_grid) rtol = rtol
    @test Array(vars_h.dynamics.div_mean) ≈ Array(vars_σ.dynamics.div_mean) rtol = rtol
    @test Array(vars_h.dynamics.div_sum_above) ≈ Array(vars_σ.dynamics.div_sum_above) rtol = rtol
    @test Array(vars_h.dynamics.pres_flux_sum_above) ≈ Array(vars_σ.dynamics.pres_flux_sum_above) rtol = rtol

    # the hybrid correction Ĉ = Σ_k (δ_k - Δσ_k) D_k is identically zero when transition ≡ 1
    # (ΔA_k ≡ 0), exactly (not just approximately)
    @test all(Array(vars_h.dynamics.div_mean_correction) .== 0)

    # vertical_velocity!
    @test Array(vars_h.dynamics.w) ≈ Array(vars_σ.dynamics.w) rtol = rtol

    # item 5 (partial, see also the standalone testset below): w is zero at the bottom
    # interface (never explicitly advected across the surface) for both coordinate types
    @test all(Array(vars_σ.dynamics.w)[:, end] .== 0)
    @test all(Array(vars_h.dynamics.w)[:, end] .== 0)

    # vertical_advection!, temperature_grid_tendency!: final grid tendencies of u, v,
    # temperature (vertical_advection! writes into these with -=, temperature_grid_tendency!
    # adds the adiabatic conversion and horizontal advection terms on top)
    for name in (:u, :v, :temperature)
        tend_σ = SpeedyWeather.get_tendency_step(getfield(vars_σ.tendencies.grid, name), sigma_model.time_stepping, SpeedyWeather.DynamicalCore())
        tend_h = SpeedyWeather.get_tendency_step(getfield(vars_h.tendencies.grid, name), hybrid_model.time_stepping, SpeedyWeather.DynamicalCore())
        @test Array(tend_h) ≈ Array(tend_σ) rtol = rtol
    end

    # surface_pressure_grid_tendency!
    p_tend_σ = SpeedyWeather.get_tendency_step(vars_σ.tendencies.grid.pressure, sigma_model.time_stepping, SpeedyWeather.DynamicalCore())
    p_tend_h = SpeedyWeather.get_tendency_step(vars_h.tendencies.grid.pressure, hybrid_model.time_stepping, SpeedyWeather.DynamicalCore())
    @test Array(p_tend_h) ≈ Array(p_tend_σ) rtol = rtol
end

@testset "mass consistency of the layer weights" begin
    # item 5: Σ_k δ_k(pₛ) ≡ Σ_k Δp_k/pₛ = 1 for any pₛ (mass conservation of the coordinate
    # definition itself), for both SigmaCoordinates and a genuine hybrid coordinate.
    nlayers = 8
    spectral_grid = SpectralGrid(nlayers = nlayers)
    σ_half = SpeedyWeather.sigma_half_spacing(nlayers)
    σ = SigmaCoordinates(spectral_grid, σ_half)
    H = CubicSigmaPressureCoordinates(spectral_grid, σ_half; pressure_only_above = 0.2, σ_only_below = 0.8)

    for coord in (σ, H), pₛ in (0.7e5, 1.0e5, 1.3e5)
        total = sum(SpeedyWeather.pressure_thickness_ratio(k, pₛ, coord) for k in 1:nlayers)
        @test total ≈ 1 rtol = 1.0e-5
    end
end

@testset "PrimitiveDryModel + CubicSigmaPressureCoordinates + HeldSuarez stays finite" begin
    # item 6: the case that would have caught the NaN bug -- a genuine (non-degenerate)
    # pressure/sigma transition, not the pure-sigma transition = _ -> 1 used everywhere else
    # before this PR.
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    coord = CubicSigmaPressureCoordinates(spectral_grid; pressure_only_above = 0.2, σ_only_below = 0.8)
    model = PrimitiveDryModel(
        spectral_grid;
        geometry = Geometry(spectral_grid; vertical_coordinates = coord),
        orography = EarthOrography(spectral_grid),
        forcing = HeldSuarez(spectral_grid),
        drag = LinearDrag(spectral_grid),
        dynamics_only = true,
    )
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation; period = Day(1))

    @test model.feedback.nans_detected == false
    @test all(isfinite, Array(simulation.variables.prognostic.vorticity))
    @test all(isfinite, Array(simulation.variables.prognostic.divergence))
    @test all(isfinite, Array(simulation.variables.prognostic.temperature))
    @test all(isfinite, Array(simulation.variables.prognostic.pressure))
end
