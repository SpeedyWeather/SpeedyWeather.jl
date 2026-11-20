using Random

# Implicit vertical diffusion of a single column on plain 1×nlayers arrays, returns the tendency
function diffuse_column(
        x, K, density_over_distance, coordinates;
        surface_pressure = 1.0e5, Δt = 1800.0, gravity = 9.81, tendency = zero(x),
        geopotential_term = nothing, couple = false,
    )
    row(v) = reshape(collect(v), 1, :)
    tend = row(tendency)
    upper, rhs = zero(tend), zero(tend)
    SpeedyWeather.implicit_vertical_diffusion!(
        1, tend, row(x), geopotential_term, row(K), row(density_over_distance),
        upper, rhs, surface_pressure, Δt, gravity, coordinates, couple,
    )
    return vec(tend)
end

# explicit diffusion operator, flux form with no flux through top and surface
function explicit_diffusion(x, K, density_over_distance, coordinates; surface_pressure = 1.0e5, gravity = 9.81)
    nlayers = length(x)
    a = [k < nlayers ? density_over_distance[k] * K[k] : 0.0 for k in 1:nlayers]   # conductance below k
    flux_below = [k < nlayers ? a[k] * (x[k] - x[k + 1]) : 0.0 for k in 1:nlayers]  # downward flux at k+½
    flux_above = [k > 1 ? flux_below[k - 1] : 0.0 for k in 1:nlayers]
    Δp = [SpeedyWeather.pressure_thickness(k, surface_pressure, coordinates) for k in 1:nlayers]
    return gravity * (flux_above - flux_below) ./ Δp
end

# ρ/Δz on the interfaces of an isothermal column at temperature T on the given coordinates
function isothermal_density_over_distance(coordinates, nlayers; surface_pressure = 1.0e5, T = 260.0, R = 287.0, gravity = 9.81)
    p_full = [SpeedyWeather.pressure(k, surface_pressure, coordinates) for k in 1:nlayers]
    return map(1:nlayers) do k
        k == nlayers && return 0.0
        p_half = SpeedyWeather.pressure_half(k + 1, surface_pressure, coordinates)
        Δz = R * T / gravity * log(p_full[k + 1] / p_full[k])
        p_half / (R * T * Δz)
    end
end

@testset "Implicit vertical diffusion operator" begin
    nlayers = 24
    spectral_grid = SpectralGrid(NF = Float64, nlayers = nlayers)
    pₛ = 0.9e5
    coordinate_systems = (
        FriersonSigmaCoordinates(spectral_grid),
        SigmaPressureCoordinates(spectral_grid),
    )

    @testset "Conservation (T1) for $(nameof(typeof(coordinates)))" for coordinates in coordinate_systems
        Δp = [SpeedyWeather.pressure_thickness(k, pₛ, coordinates) for k in 1:nlayers]
        rng = Xoshiro(123)
        for couple in (false, true), Δt in (0.0, 300.0, 1.0e5)
            x = 280 .+ 10 * randn(rng, nlayers)
            K = 100 * rand(rng, nlayers)
            ρ_Δz = 1.0e-2 * rand(rng, nlayers)
            source = 1.0e-4 * randn(rng, nlayers)
            tend = diffuse_column(x, K, ρ_Δz, coordinates; surface_pressure = pₛ, Δt, tendency = source, couple)
            @test sum(tend .* Δp) ≈ sum(source .* Δp) atol = 1.0e-12 * sum(abs.(x .* Δp))
        end
    end

    @testset "Well-mixed and adiabatic columns (T2)" begin
        coordinates = first(coordinate_systems)
        K = fill(50.0, nlayers)
        ρ_Δz = isothermal_density_over_distance(coordinates, nlayers)

        # constant column: zero tendency
        @test all(diffuse_column(fill(3.0, nlayers), K, ρ_Δz, coordinates) .== 0)

        # adiabatic column, dry static energy cₚT + Φ constant: zero temperature tendency
        cₚ = 1004.0
        Φ = range(3.0e5, 400.0, length = nlayers)
        T = 300 .- Φ / cₚ
        tend = diffuse_column(T, K, ρ_Δz, coordinates; geopotential_term = (reshape(collect(Φ), 1, :), 1 / cₚ))
        @test all(abs.(tend) .< 1.0e-12)
    end

    @testset "Units: decay of a cosine mode (T3)" begin
        # incompressible column: equally spaced σ, constant ρ and Δz, so L = K ∂²/∂z²
        coordinates = SigmaCoordinates(spectral_grid; σ_half = collect(range(0, 1, length = nlayers + 1)))
        g = 9.81
        Δz = 100.0
        H = nlayers * Δz
        Δp = SpeedyWeather.pressure_thickness(1, pₛ, coordinates)
        ρ = Δp / (g * Δz)
        K = 20.0
        ρ_Δz = [k < nlayers ? ρ / Δz : 0.0 for k in 1:nlayers]
        x = [cos(π * (k - 1 / 2) / nlayers) for k in 1:nlayers]

        tend = diffuse_column(x, fill(K, nlayers), ρ_Δz, coordinates; surface_pressure = pₛ, Δt = 0.0, gravity = g)
        decay_rate_discrete = 4K / Δz^2 * sin(π / (2nlayers))^2
        @test tend ≈ -decay_rate_discrete * x rtol = 1.0e-10
        @test tend ≈ -K * (π / H)^2 * x rtol = 1.0e-2
    end

    @testset "Unconditional stability and monotonicity (T4)" begin
        coordinates = first(coordinate_systems)
        ρ_Δz = isothermal_density_over_distance(coordinates, nlayers)
        rng = Xoshiro(1)
        x = 280 .+ 10 * randn(rng, nlayers)
        Δt = 3600.0
        tend = diffuse_column(x, fill(1.0e4, nlayers), ρ_Δz, coordinates; Δt)
        x_new = x + Δt * tend
        @test all(isfinite, x_new)
        @test minimum(x) <= minimum(x_new) && maximum(x_new) <= maximum(x)  # no new extrema

        # mass-weighted variance decreases (diffusion length √(KΔt) ≈ 6 km, so not fully mixed)
        Δp = [SpeedyWeather.pressure_thickness(k, 1.0e5, coordinates) for k in 1:nlayers]
        variance(x) = sum(Δp .* (x .- sum(Δp .* x) / sum(Δp)) .^ 2) / sum(Δp)
        @test variance(x_new) < variance(x) / 2
    end

    @testset "Implicit converges to explicit for small Δt (T5)" begin
        coordinates = first(coordinate_systems)
        ρ_Δz = isothermal_density_over_distance(coordinates, nlayers)
        rng = Xoshiro(2)
        x = 280 .+ 10 * randn(rng, nlayers)
        K = 30 * rand(rng, nlayers)
        explicit = explicit_diffusion(x, K, ρ_Δz, coordinates)
        @test diffuse_column(x, K, ρ_Δz, coordinates; Δt = 0.0) ≈ explicit rtol = 1.0e-12
        error_1 = maximum(abs, diffuse_column(x, K, ρ_Δz, coordinates; Δt = 1.0) - explicit)
        error_10 = maximum(abs, diffuse_column(x, K, ρ_Δz, coordinates; Δt = 10.0) - explicit)
        @test error_1 < 1.0e-2 * maximum(abs, explicit)
        @test error_10 / error_1 ≈ 10 rtol = 0.1      # first order in Δt
    end
end

# model with vertical diffusion after a day, so that low-level winds and a boundary layer develop
function vertical_diffusion_simulation(Model = PrimitiveWetModel; nlayers = 8, kwargs...)
    spectral_grid = SpectralGrid(truncation = 32; nlayers)
    model = Model(spectral_grid; vertical_diffusion = BulkRichardsonDiffusion(spectral_grid; kwargs...))
    simulation = initialize!(model)
    run!(simulation, period = Day(1))
    return simulation
end

# recompute the diffusivities on all grid points
function vertical_diffusivity!(simulation)
    (; model, variables) = simulation
    for ij in 1:model.geometry.npoints
        SpeedyWeather.vertical_diffusivity!(ij, variables, model.vertical_diffusion, model)
    end
    p = variables.parameterizations
    return copy(p.vertical_diffusivity_momentum.data), copy(p.vertical_diffusivity_heat.data), copy(p.boundary_layer_height.data)
end

# height [m] of the interfaces k+½ above the surface
function interface_heights(simulation)
    (; model, variables) = simulation
    Φ = variables.dynamics.geopotential.data
    z = (Φ[:, 1:(end - 1)] + Φ[:, 2:end]) / 2 / model.planet.gravity .- model.orography.orography.data
    return hcat(z, zeros(size(z, 1)))
end

@testset "BulkRichardsonDiffusion diffusivity" begin
    simulation = vertical_diffusion_simulation()
    (; model, variables) = simulation
    diffusion = model.vertical_diffusion
    (; background_diffusivity, background_height, max_diffusivity) = diffusion
    nlayers = model.geometry.nlayers
    g = model.planet.gravity

    @testset "Bounds" begin
        K_momentum, K_heat, h = vertical_diffusivity!(simulation)
        z = interface_heights(simulation)
        K_background = background_diffusivity * exp.(-z / background_height)
        for K in (K_momentum, K_heat)
            @test all(isfinite, K)
            @test all(K[:, end] .== 0)                                      # no flux through the surface
            @test all(K[:, 1:(end - 1)] .>= K_background[:, 1:(end - 1)] * (1 - 1.0e-5))
            @test all(K .<= max_diffusivity)
        end
        @test all(K_momentum .>= K_heat)       # longer tail for momentum

        # boundary layer height between lowermost and uppermost layer
        z_layers = variables.dynamics.geopotential.data / g .- model.orography.orography.data
        @test all(z_layers[:, end] .<= h .<= z_layers[:, 1])
    end

    @testset "Response to the surface temperature (T6)" begin
        air_temperature = SpeedyWeather.get_prognostic_step(variables.grid.temperature, model.time_stepping, diffusion).data[:, end]
        drag = model.boundary_layer.drag

        function diffusivity(ΔT, surface_wind_speed)
            variables.prognostic.land.soil_temperature.data .= air_temperature .+ ΔT
            variables.prognostic.ocean.sea_surface_temperature.data .= air_temperature .+ ΔT
            variables.parameterizations.surface_wind_speed.data .= surface_wind_speed
            for ij in eachindex(air_temperature)
                SpeedyWeather.boundary_layer_drag!(ij, variables, drag, model)
            end
            return vertical_diffusivity!(simulation)
        end

        K_warm, _, h = diffusivity(20, 10)              # unstable surface layer, Ri < 0
        K_cold, _, h_cold = diffusivity(-10, 10)        # stable surface layer, 0 < Ri < Ri_c
        K_very_cold, _, _ = diffusivity(-20, 0.5)       # Ri ≥ Ri_c at the surface, no boundary layer profile

        @test h == h_cold                               # h from the atmospheric column only
        @test all(K_warm .>= K_cold .>= K_very_cold)

        # the profile only acts where the lowermost interface is within the boundary layer
        z = interface_heights(simulation)
        inside = h .> z[:, end - 1]
        @test count(inside) > 0
        @test sum(K_warm[inside, end - 1]) > sum(K_cold[inside, end - 1]) > sum(K_very_cold[inside, end - 1])

        # without the profile only the local and background diffusivities remain
        @test all(K_very_cold[:, end - 1] .>= background_diffusivity * exp.(-z[:, end - 1] / background_height) * (1 - 1.0e-5))
    end

    @testset "Independence of the surface height (T7)" begin
        K_momentum, K_heat, h = vertical_diffusivity!(simulation)
        Φ = variables.dynamics.geopotential.data
        orography = model.orography.orography.data
        Φ .+= g * 4000
        orography .+= 4000
        K_momentum_4km, K_heat_4km, h_4km = vertical_diffusivity!(simulation)
        Φ .-= g * 4000
        orography .-= 4000
        @test K_momentum_4km ≈ K_momentum rtol = 1.0e-3
        @test K_heat_4km ≈ K_heat rtol = 1.0e-3
        @test h_4km ≈ h rtol = 1.0e-3
    end
end

@testset "BulkRichardsonDiffusion mixing" begin
    simulation = vertical_diffusion_simulation()
    (; model, variables) = simulation
    diffusion = model.vertical_diffusion
    nlayers = model.geometry.nlayers
    npoints = model.geometry.npoints
    (; vertical_coordinates) = model.geometry
    g = model.planet.gravity
    cₚ = model.atmosphere.heat_capacity
    Δt = SpeedyWeather.default_time_step(model.time_stepping)
    T = SpeedyWeather.get_prognostic_step(variables.grid.temperature, model.time_stepping, diffusion).data
    Φ = variables.dynamics.geopotential.data
    pₛ = variables.parameterizations.surface_pressure.data

    # temperature tendency of the vertical diffusion of dry static energy for all columns
    function temperature_tendency(; source = zeros(npoints, nlayers), couple = true)
        _, K_heat, _ = vertical_diffusivity!(simulation)
        ρ_Δz = variables.scratch.grid.density_over_distance.data
        return reduce(
            vcat, map(1:npoints) do ij
                diffuse_column(
                    T[ij, :], K_heat[ij, :], ρ_Δz[ij, :], vertical_coordinates;
                    surface_pressure = pₛ[ij], Δt, gravity = g, tendency = source[ij, :], couple,
                    geopotential_term = (Φ[ij:ij, :], 1 / cₚ),
                )'
            end
        )
    end

    @testset "Superadiabatic layer is mixed (T8)" begin
        k = nlayers ÷ 2         # mid-troposphere
        T_original = copy(T)
        T[:, k] .-= 30          # superadiabatic interface k+½: dry static energy decreasing upwards
        Δs(T) = cₚ * (T[:, k] - T[:, k + 1]) + Φ[:, k] - Φ[:, k + 1]
        Δs_initial = Δs(T)
        @test all(Δs_initial .< 0)

        _, K_heat, _ = vertical_diffusivity!(simulation)
        @test all(K_heat[:, k] .> 10)      # [m²/s] unstable: local diffusivity also above the boundary layer

        # a few time steps with Φ held fixed, adjusts within hours, monotonically
        Δs_previous = Δs_initial
        for _ in 1:3
            T .+= Δt * temperature_tendency()
            @test all(Δs_previous .< Δs(T) .< 0)
            Δs_previous = Δs(T)
        end
        @test all(Δs(T) .> Δs_initial / 2)      # at least halfway towards neutral
        T .= T_original
    end

    @testset "Surface flux is spread over the boundary layer (T9)" begin
        # heating of the lowermost layer as by a surface sensible heat flux
        source = zeros(npoints, nlayers)
        source[:, end] .= 1.0e-4            # [K/s]
        coupled = temperature_tendency(; source, couple = true)
        uncoupled = temperature_tendency(; source, couple = false)
        spread = coupled - uncoupled

        Δp = [SpeedyWeather.pressure_thickness(k, pₛ[ij], vertical_coordinates) for ij in 1:npoints, k in 1:nlayers]
        @test all(spread[:, end] .< 0)                      # less heating of the lowermost layer
        @test all(spread[:, end - 1] .> 0)                  # more heating of the layer above
        @test all(abs.(sum(spread .* Δp, dims = 2)) .< 1.0e-6 * 1.0e-4 * Δp[:, end])  # same column heating
    end
end

@testset "BulkRichardsonDiffusion in a model" begin
    @testset for (Model, nlayers) in ((PrimitiveWetModel, 8), (PrimitiveWetModel, 24), (PrimitiveDryModel, 24))
        spectral_grid = SpectralGrid(truncation = 32; nlayers)
        model = Model(spectral_grid)
        @test model.vertical_diffusion isa BulkRichardsonDiffusion      # default
        simulation = initialize!(model)
        run!(simulation, period = Day(2))
        @test !simulation.model.feedback.nans_detected
        (; vertical_diffusivity_momentum, vertical_diffusivity_heat, boundary_layer_height) = simulation.variables.parameterizations
        @test all(isfinite, vertical_diffusivity_momentum)
        @test all(isfinite, vertical_diffusivity_heat)
        @test all(boundary_layer_height .> 0)
    end

    # disabled diffusion leaves the tendencies unchanged
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    diffusion = BulkRichardsonDiffusion(spectral_grid; diffuse_momentum = false, diffuse_static_energy = false, diffuse_humidity = false)
    model = PrimitiveWetModel(spectral_grid; vertical_diffusion = diffusion)
    simulation = initialize!(model)
    run!(simulation, steps = 3)
    @test all(simulation.variables.parameterizations.vertical_diffusivity_heat .== 0)   # never computed

    # friction velocity and surface Richardson number from the boundary layer, both zero without its
    # drag and surface wind, which switches the boundary-layer profile of the diffusivity off
    (; variables) = simulation
    ΔΦ₀ = 500.0f0       # geopotential of the lowermost layer relative to the surface
    friction_velocity, Ri_surface = SpeedyWeather.surface_layer(1, ΔΦ₀, variables, diffusion, model)
    @test friction_velocity > 0
    @test isfinite(Ri_surface)
    no_boundary_layer = (; parameterizations = (; variables.parameterizations.vertical_diffusivity_momentum))
    @test SpeedyWeather.surface_layer(1, ΔΦ₀, no_boundary_layer, diffusion, model) === (0.0f0, 0.0f0)
end
