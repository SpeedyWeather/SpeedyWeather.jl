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

    # eastward at the equator: the great-circle angle is u*Δt/radius
    radius = 6.371e6
    Δt_r = Float32(SpeedyWeather.trajectory_time_step(3600.0, radius))
    u = 10.0f0      # m/s eastward for an hour at the equator

    lond, latd = displace(0.0f0, 0.0f0, u, 0.0f0, Δt_r)
    @test latd ≈ 0 atol = 1.0f-5
    @test lond ≈ Float32(u * 3600 / radius * 180 / π) rtol = 1.0f-4

    # longitudes wrap into [0, 360)
    lond, latd = displace(359.0f0, 0.0f0, 1000.0f0, 0.0f0, Δt_r)
    @test 0 <= lond < 360

    # NO POLE SINGULARITY. The lat/lon form needs dlon = u*Δt/(radius*cos(lat)), which blows up
    # at high latitude: at 89.28˚N (T128 outermost ring) cos(lat) = 0.0125, so a 50 m/s wind over
    # 45 min would give a 97˚ longitude jump. The Cartesian form must stay on the sphere.
    Δt_45min = Float32(SpeedyWeather.trajectory_time_step(2700.0, radius))
    lond, latd = displace(0.0f0, 89.28f0, 50.0f0, 0.0f0, Δt_45min)
    @test -90 <= latd <= 90
    @test 0 <= lond < 360
    # a purely zonal wind at 89.28˚N sweeps a great circle, which *must* carry it over the pole
    # rather than around a latitude circle; either way the displacement is bounded by u*Δt/radius
    angular_distance = acosd(
        clamp(
            sind(89.28f0) * sind(latd) + cosd(89.28f0) * cosd(latd) * cosd(lond - 0.0f0),
            -1.0f0, 1.0f0
        )
    )
    @test angular_distance ≈ rad2deg(50 * 2700 / radius) rtol = 1.0f-3

    # Exact round trip along the equator, where a zonal wind stays zonal and the local east
    # vector is parallel-transported, so (u, v) mean the same thing at both ends. Away from the
    # equator the local frame rotates along the path, so a round trip with the *same* (u, v)
    # components is not expected to close exactly — that is geometry, not an error.
    lon1, lat1 = displace(35.0f0, 0.0f0, 30.0f0, 0.0f0, Δt_45min)
    lon0, lat0 = displace(lon1, lat1, 30.0f0, 0.0f0, -Δt_45min)
    @test lon0 ≈ 35.0f0 atol = 1.0f-4
    @test lat0 ≈ 0.0f0 atol = 1.0f-4
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
    @test haskey(vars.dynamics, :semi_lagrangian)

    # departure points are grid fields, one per arrival grid cell
    @test vars.dynamics.semi_lagrangian.departure_lond isa RingGrids.AbstractField

    # the time stepper must NOT reach into the diffusion component: that choice belongs to
    # HyperDiffusion alone, and SemiLagrangian works with either damping factor
    @test model.horizontal_diffusion.exponential == false
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
    (; departure_lond, departure_latd) = vars.dynamics.semi_lagrangian
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

@testset "SemiLagrangian: previous grid step initialized" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    model = BarotropicModel(
        spectral_grid,
        time_stepping = SemiLagrangian(spectral_grid),
        forcing = nothing, drag = nothing, random_process = nothing,
    )
    simulation = initialize!(model)
    set!(simulation, vorticity = (lon, lat, σ) -> 1.0e-5 * exp(-((lon - 180)^2 + lat^2) / 200))
    SpeedyWeather.transform!(simulation.variables, model, initialize = true)

    # step 1 (previous) must be a copy of step 2 (current) after initialization, otherwise the
    # SETTLS extrapolation 3/2 uⁿ - 1/2 uⁿ⁻¹ would read zeros on the first step
    u = simulation.variables.grid.u
    @test view(u.data, :, 1, 1) == view(u.data, :, 1, 2)
    @test maximum(abs, view(u.data, :, 1, 1)) > 0

    # Leapfrog on a 2D model keeps one grid step, the hook must be a no-op there
    model_lf = BarotropicModel(spectral_grid)
    sim_lf = initialize!(model_lf)
    @test size(sim_lf.variables.grid.u, 3) == 1
    @test SpeedyWeather.transform!(sim_lf.variables, model_lf, initialize = true) === nothing
end

@testset "SemiLagrangian: displace compiles away" begin
    # the Particle round trip in `displace` must not allocate, it runs per grid point per iteration
    @test Base.return_types(SpeedyWeather.displace, NTuple{5, Float32}) == [Tuple{Float32, Float32}]
end

@testset "SemiLagrangian: interpolator is selectable" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)

    # cubic is the default, since bilinear-class interpolation damps the transported field
    # every time step (see RingGrids/test/interpolation_cubic.jl)
    default_stepper = SemiLagrangian(spectral_grid)
    @test default_stepper.interpolator isa RingGrids.CubicInterpolator

    # the trajectory only has to place the departure point, so it defaults to the cheap one:
    # the damping that matters accumulates on the transported field, not on the winds
    @test default_stepper.trajectory_interpolator isa RingGrids.AnvilInterpolator

    for Interpolator in (RingGrids.CubicInterpolator, RingGrids.AnvilInterpolator)
        time_stepping = SemiLagrangian(spectral_grid; Interpolator)
        @test time_stepping.interpolator isa Interpolator

        model = BarotropicModel(
            spectral_grid; time_stepping,
            forcing = nothing, drag = nothing, random_process = nothing,
        )
        simulation = initialize!(model)
        set!(simulation, vorticity = (lon, lat, σ) -> 1.0e-5 * exp(-((lon - 180)^2 + lat^2) / 200))
        run!(simulation, period = Hour(6))
        @test all(isfinite, simulation.variables.grid.vorticity)
    end
end

@testset "SemiLagrangian: trajectory interpolator and iteration count" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    C, A = RingGrids.CubicInterpolator, RingGrids.AnvilInterpolator

    function run6h(; kwargs...)
        time_stepping = SemiLagrangian(spectral_grid; kwargs...)
        model = BarotropicModel(
            spectral_grid; time_stepping,
            forcing = nothing, drag = nothing, random_process = nothing,
        )
        simulation = initialize!(model)
        set!(simulation, vorticity = (lon, lat, σ) -> 1.0e-5 * exp(-((lon - 180)^2 + lat^2) / 200))
        run!(simulation, period = Hour(6))
        return simulation.variables
    end

    # both interpolators are independently selectable
    for Interpolator in (C, A), TrajectoryInterpolator in (C, A)
        time_stepping = SemiLagrangian(spectral_grid; Interpolator, TrajectoryInterpolator)
        @test time_stepping.interpolator isa Interpolator
        @test time_stepping.trajectory_interpolator isa TrajectoryInterpolator
    end

    # a cheap trajectory interpolator must not change the answer much: it only places the
    # departure point, and the winds are rebuilt from the spectral state every step anyway
    vars_cheap = run6h(Interpolator = C, TrajectoryInterpolator = A)
    vars_full = run6h(Interpolator = C, TrajectoryInterpolator = C)
    @test maximum(abs, vars_cheap.grid.vorticity .- vars_full.grid.vorticity) <
        0.05 * maximum(abs, vars_full.grid.vorticity)

    # the trajectory loop skips the redundant first interpolation (departure = arrival, where the
    # wind is already known), so n_iterations = 1 is a valid Euler-backward trajectory
    for n_iterations in (1, 2, 4)
        vars = run6h(; n_iterations)
        @test all(isfinite, vars.grid.vorticity)
        (; departure_lond, departure_latd) = vars.dynamics.semi_lagrangian
        @test all(lon -> 0 <= lon < 360, departure_lond)
        @test all(lat -> -90 <= lat <= 90, departure_latd)
    end
end
