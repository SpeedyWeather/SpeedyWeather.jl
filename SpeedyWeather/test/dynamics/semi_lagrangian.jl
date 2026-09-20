@testset "φ₁ exponential integrator function" begin
    φ₁ = SpeedyWeather.φ₁

    # φ₁(0) = 1 exactly, no 0/0
    @test φ₁(0.0f0) == 1
    @test φ₁(0.0) == 1

    # the defining property: 1 + z*φ₁(z) == exp(z), which is what makes the
    # (tendency + expl*var)*impl step exact for the diffusion operator
    for NF in (Float32, Float64)
        for z in NF.((-1.0e-9, -1.0e-4, -0.125, -1, -5))
            @test 1 + z * φ₁(z) ≈ exp(z) rtol = 10eps(NF)
        end
    end

    # and it stays close to the backward Euler factor for small |z|, where both → 1
    @test φ₁(-1.0e-6) ≈ 1 / (1 + 1.0e-6) rtol = 1.0e-5
end

@testset "HyperDiffusion exponential option" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)

    # default is backward Euler, opt-in is the exponential
    model_euler = BarotropicModel(spectral_grid)
    model_exp = BarotropicModel(spectral_grid, horizontal_diffusion = HyperDiffusion(spectral_grid, exponential = true))
    initialize!(model_euler)
    initialize!(model_exp)

    expl = model_euler.horizontal_diffusion.expl
    Δt = SpeedyWeather.default_time_step(model_euler.time_stepping) / model_euler.planet.radius

    # the explicit part is untouched by the option, only `impl` differs
    @test model_euler.horizontal_diffusion.expl == model_exp.horizontal_diffusion.expl

    for l in 1:(spectral_grid.truncation - 1)
        z = Δt * expl[l, 1]
        @test model_euler.horizontal_diffusion.impl[l, 1] ≈ 1 / (1 - z)
        @test model_exp.horizontal_diffusion.impl[l, 1] ≈ SpeedyWeather.φ₁(z)
    end
end

@testset "SemiLagrangian: displacement on the sphere" begin
    displace = SpeedyWeather.displace

    # zero velocity leaves a point where it is
    @test all(displace(120.0f0, 30.0f0, 0.0f0, 0.0f0, 1.0f0) .≈ (120.0f0, 30.0f0))

    # Δt_deg converts [m/s] to [˚]: u*Δt/radius*180/π
    radius = 6.371e6
    Δt_deg = Float32(SpeedyWeather.trajectory_time_step(3600.0, radius))
    u = 10.0f0      # m/s eastward for an hour at the equator

    lond, latd = displace(0.0f0, 0.0f0, u, 0.0f0, Δt_deg)
    @test latd ≈ 0 atol = 1.0f-5
    @test lond ≈ Float32(u * 3600 / radius * 180 / π) rtol = 1.0f-4

    # longitudes wrap into [0, 360)
    lond, latd = displace(359.0f0, 0.0f0, 1000.0f0, 0.0f0, Δt_deg)
    @test 0 <= lond < 360
end

@testset "SemiLagrangian: step structure" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    time_stepping = SemiLagrangian(spectral_grid)

    # two-time-level: one spectral state, two grid levels for the wind extrapolation
    @test SpeedyWeather.prognostic_spectral_steps(time_stepping) == 1
    @test SpeedyWeather.prognostic_grid_steps(time_stepping, BarotropicModel(spectral_grid)) == 2

    # the trajectory shift replaces the vorticity flux
    @test SpeedyWeather.advection_factor(time_stepping) == 0
    @test SpeedyWeather.advection_factor(Leapfrog(spectral_grid)) == 1

    model = BarotropicModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    vars = simulation.variables

    @test size(vars.prognostic.vorticity, 3) == 1   # spectral: 1 step
    @test size(vars.grid.u, 3) == 2                 # grid: 2 steps
    @test haskey(vars.scratch, :semi_lagrangian)

    # opted into the exponential diffusion
    @test model.horizontal_diffusion.exponential
end

@testset "SemiLagrangian: barotropic model runs" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    model = BarotropicModel(
        spectral_grid,
        time_stepping = SemiLagrangian(spectral_grid),
        forcing = nothing, drag = nothing, random_process = nothing,
    )
    simulation = initialize!(model)
    set!(simulation, vorticity = (lon, lat, σ) -> 1.0e-5 * exp(-((lon - 180)^2 + lat^2) / 200))

    run!(simulation, period = Hour(6))
    vars = simulation.variables

    @test all(isfinite, vars.grid.vorticity)
    @test all(isfinite, vars.grid.u)

    # departure points stay on the sphere
    (; departure_lond, departure_latd) = vars.scratch.semi_lagrangian
    @test all(lon -> 0 <= lon < 360, departure_lond)
    @test all(lat -> -90 <= lat <= 90, departure_latd)

    # without forcing or drag the flow should not spin up out of nowhere
    @test maximum(abs, vars.grid.vorticity) < 1.0e-4
end

@testset "SemiLagrangian: only barotropic supported" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 4)
    model = PrimitiveDryModel(spectral_grid, time_stepping = SemiLagrangian(spectral_grid))
    @test_throws ArgumentError initialize!(model)
end
