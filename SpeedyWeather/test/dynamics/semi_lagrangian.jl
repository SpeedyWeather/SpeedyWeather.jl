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
    @test default_stepper.Locator === RingGrids.CubicLocator

    # the trajectory only has to place the departure point, so it defaults to the cheap one:
    # the damping that matters accumulates on the transported field, not on the winds
    @test default_stepper.TrajectoryLocator === RingGrids.AnvilLocator

    # the geometry is constant and stays in the model, the locators vary and live in Variables
    @test default_stepper.geometry isa RingGrids.GridGeometry
    @test !any(
        f -> fieldtype(typeof(default_stepper), f) <: RingGrids.AbstractLocator,
        fieldnames(typeof(default_stepper))
    )

    for Interpolator in (RingGrids.CubicInterpolator, RingGrids.AnvilInterpolator)
        time_stepping = SemiLagrangian(spectral_grid; Interpolator)
        @test time_stepping.Locator === RingGrids.Locator(Interpolator)

        model = BarotropicModel(
            spectral_grid; time_stepping,
            forcing = nothing, drag = nothing, random_process = nothing,
        )
        simulation = initialize!(model)
        @test simulation.variables.dynamics.semi_lagrangian.locator isa RingGrids.Locator(Interpolator)
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
        model = BarotropicModel(
            spectral_grid; time_stepping,
            forcing = nothing, drag = nothing, random_process = nothing,
        )
        sl = initialize!(model).variables.dynamics.semi_lagrangian
        @test sl.locator isa RingGrids.Locator(Interpolator)
        @test sl.trajectory_locator isa RingGrids.Locator(TrajectoryInterpolator)
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

@testset "SemiLagrangian: trajectory winds combine in a common frame" begin
    # `u, v` are components in the *local* east/north frame. Averaging the arrival and departure
    # winds componentwise is only valid where the two frames coincide; the mismatch is the
    # meridian convergence and is O(1) next to the poles, where it made the fixed-point iteration
    # diverge instead of converge. Both winds are therefore mapped to 3D Cartesian first.
    cv = SpeedyWeather.cartesian_velocity

    # a purely zonal wind at two points half a polar ring apart points in nearly opposite
    # directions in 3D, even though both have the same (u, v) = (50, 0)
    Va = cv(0.0, 89.0, 50.0, 0.0)
    Vd = cv(180.0, 89.0, 50.0, 0.0)
    @test Va[1] ≈ -Vd[1] atol = 1.0e-10
    @test Va[2] ≈ -Vd[2] atol = 1.0e-10
    # componentwise the mean would be (50, 0), in the common frame it is ~zero
    @test maximum(abs, (Va .+ Vd) ./ 2) < 1.0e-10

    # the Cartesian velocity is tangent to the sphere: V ⟂ r
    for (lond, latd, u, v) in ((0.0, 0.0, 30.0, -20.0), (73.0, 89.5, -40.0, 10.0), (250.0, -60.0, 5.0, 5.0))
        r = (cosd(latd) * cosd(lond), cosd(latd) * sind(lond), sind(latd))
        V = cv(lond, latd, u, v)
        @test abs(sum(V .* r)) < 1.0e-12
        @test sqrt(sum(abs2, V)) ≈ sqrt(u^2 + v^2)        # and norm-preserving
    end
end

@testset "SemiLagrangian: trajectory covers the full time step" begin
    # x_d = x_a - Δt (V(x_a) + V(x_d))/2: the mean wind already carries the 1/2, so the
    # displacement must be the full Δt. Halving it once made the scheme advect at half speed.
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    time_stepping = SemiLagrangian(spectral_grid)
    model = BarotropicModel(
        spectral_grid; time_stepping,
        forcing = nothing, drag = nothing, random_process = nothing,
    )
    simulation = initialize!(model)
    vars = simulation.variables

    # a uniform zonal wind: every trapezoidal iterate is exact, so the departure point is
    # exactly one `displace` of -Δt away from the arrival point
    U = 50.0f0
    for step in (1, 2)
        fill!(SpeedyWeather.field_view(vars.grid.u, :, 1, step).data, U)
        fill!(SpeedyWeather.field_view(vars.grid.v, :, 1, step).data, 0)
    end
    SpeedyWeather.departure_points!(vars, time_stepping, model)

    Δt_r = time_stepping.Δt_trajectory[]
    sl = vars.dynamics.semi_lagrangian
    londs, latds = model.geometry.londs, model.geometry.latds

    # Compare in 3D to sidestep the longitude wrap. Restricted to |lat| <= 45˚: a uniform (u, v)
    # is not a great-circle flow, so arrival and departure frames genuinely disagree, and the
    # reference below is only exact where they coincide. The residual grows smoothly with
    # latitude (9e-5 at 15˚, 2.7e-4 at 45˚, 7e-3 at the pole).
    p(lond, latd) = (cosd(latd) * cosd(lond), cosd(latd) * sind(lond), sind(latd))
    worst, halved = 0.0, 0.0
    for ij in eachindex(londs)
        abs(latds[ij]) <= 45 || continue
        arrived = p(sl.departure_lond[ij], sl.departure_latd[ij])
        worst = max(worst, maximum(abs, p(SpeedyWeather.displace(londs[ij], latds[ij], U, 0.0f0, -Δt_r)...) .- arrived))
        halved = max(halved, maximum(abs, p(SpeedyWeather.displace(londs[ij], latds[ij], U, 0.0f0, -Δt_r / 2)...) .- arrived))
    end
    @test worst < 1.0e-3

    # and the tolerance is tight enough to catch a trajectory of the wrong length
    @test halved > 20 * worst
end

@testset "SemiLagrangian: previous wind is the previous step, not the current one" begin
    # `move_prognostic_grid_variables_back!` has to run *before* the transforms overwrite the
    # current grid step. Called afterwards it copies the new wind onto the previous one, so
    # 3/2 uⁿ - 1/2 uⁿ⁻¹ silently collapses to uⁿ and the Rossby term becomes forward Euler.
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
    time_stepping = SemiLagrangian(spectral_grid)
    model = BarotropicModel(
        spectral_grid; time_stepping,
        forcing = nothing, drag = nothing, random_process = nothing,
    )
    simulation = initialize!(model)
    set!(simulation, vorticity = (lon, lat, σ) -> 1.0e-5 * exp(-((lon - 180)^2 + lat^2) / 200))
    run!(simulation, period = Hour(6))

    u = simulation.variables.grid.u
    u_old = view(u.data, :, 1, 1)
    u_new = view(u.data, :, 1, 2)
    @test u_old != u_new                                    # the two steps must differ
    @test maximum(abs, u_new .- u_old) > 1.0e-4 * maximum(abs, u_new)

    # and the extrapolated wind must actually extrapolate, i.e. leave the 1-step interval
    SpeedyWeather.extrapolate_winds!(simulation.variables, time_stepping, model)
    sl = simulation.variables.dynamics.semi_lagrangian
    @test maximum(abs, sl.u_star .- u_new) > 0
end

@testset "SemiLagrangian: Rossby waves do not amplify" begin
    # The trajectory carries the Rossby term through f(x_d) - f(x_a), so the time level of the
    # trajectory wind *is* the time discretisation of the Rossby wave. With the SETTLS
    # extrapolation this is AB2-like and near-neutral; with the wind at tⁿ it is forward Euler
    # and unconditionally unstable. Free decay of a small-amplitude m=2, n=3 mode, no forcing,
    # no drag, negligible diffusion.
    function amplification(; extrapolate_winds)
        spectral_grid = SpectralGrid(truncation = 31, nlayers = 1)
        time_stepping = SemiLagrangian(spectral_grid; extrapolate_winds)
        model = BarotropicModel(
            spectral_grid; time_stepping,
            forcing = nothing, drag = nothing, random_process = nothing,
            horizontal_diffusion = HyperDiffusion(spectral_grid, time_scale = Day(10000)),
        )
        simulation = initialize!(model)
        set!(simulation, vorticity = (λ, φ, σ) -> 1.0e-6 * cosd(φ)^2 * sind(φ) * cosd(2λ))
        e0 = Float64(sum(abs2, simulation.variables.prognostic.vorticity))
        run!(simulation, period = Day(2))
        e1 = Float64(sum(abs2, simulation.variables.prognostic.vorticity))
        return sqrt(e1 / e0)
    end

    # over 2 days the wave must keep its amplitude to within a per cent
    @test 0.99 < amplification(extrapolate_winds = true) < 1.01

    # and the un-extrapolated variant is the forward-Euler case it is documented to be: unstable
    @test amplification(extrapolate_winds = false) > 1.05
end
