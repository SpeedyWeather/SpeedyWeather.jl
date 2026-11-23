using Random

# Implicit vertical diffusion of a single column on plain 1×nlayers arrays, returns the tendency. The solver
# takes the conductance ρK/Δz of the interfaces, here from the diffusivity K and ρ/Δz
function diffuse_column(
        x, K, density_over_distance, coordinates;
        surface_pressure = 1.0e5, Δt = 1800.0, gravity = 9.81, tendency = zero(x),
        geopotential_term = nothing, couple = false,
    )
    row(v) = reshape(collect(v), 1, :)
    tend = row(tendency)
    upper, rhs = fill(NaN, size(tend)), fill(NaN, size(tend))      # scratch, to be written before read
    SpeedyWeather.implicit_vertical_diffusion!(
        1:1, ((tend, row(x), geopotential_term),), row(collect(K) .* collect(density_over_distance)),
        upper, rhs, [surface_pressure], Δt, gravity, coordinates, couple,
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

    @testset "Several variables and columns at once" begin
        coordinates = first(coordinate_systems)
        ncolumns = 7
        cₚ = 1004.0
        rng = Xoshiro(3)
        x = 280 .+ 10 * randn(rng, ncolumns, nlayers)       # with a geopotential term
        y = 1.0e-2 * rand(rng, ncolumns, nlayers)           # without
        Φ = repeat(reshape(collect(range(3.0e5, 400.0, length = nlayers)), 1, :), ncolumns) .* (1 .+ 0.1 * rand(rng, ncolumns))
        K = 100 * rand(rng, ncolumns, nlayers)
        ρ_Δz = 1.0e-2 * rand(rng, ncolumns, nlayers)
        surface_pressure = 1.0e5 .+ 1.0e3 * randn(rng, ncolumns)
        source_x = 1.0e-4 * randn(rng, ncolumns, nlayers)
        source_y = 1.0e-8 * randn(rng, ncolumns, nlayers)
        Δt, g = 1800.0, 9.81

        # both variables for the columns of a range with the solver for ranges
        function diffuse_columns(columns, couple)
            tendency_x, tendency_y = copy(source_x), copy(source_y)
            upper, rhs = fill(NaN, size(x)), fill(NaN, size(x))     # scratch, to be written before read
            fields = ((tendency_x, x, (Φ, 1 / cₚ)), (tendency_y, y, nothing))
            SpeedyWeather.implicit_vertical_diffusion!(columns, fields, K .* ρ_Δz, upper, rhs, surface_pressure, Δt, g, coordinates, couple)
            return tendency_x, tendency_y
        end

        for couple in (true, false)
            # reference: one variable and one column at a time
            reference_x = reduce(
                vcat, [
                    diffuse_column(
                        x[ij, :], K[ij, :], ρ_Δz[ij, :], coordinates; surface_pressure = surface_pressure[ij], Δt, gravity = g,
                        tendency = source_x[ij, :], geopotential_term = (Φ[ij:ij, :], 1 / cₚ), couple
                    )' for ij in 1:ncolumns
                ]
            )
            reference_y = reduce(
                vcat, [
                    diffuse_column(
                        y[ij, :], K[ij, :], ρ_Δz[ij, :], coordinates; surface_pressure = surface_pressure[ij], Δt, gravity = g,
                        tendency = source_y[ij, :], couple
                    )' for ij in 1:ncolumns
                ]
            )
            tendency_x, tendency_y = diffuse_columns(1:ncolumns, couple)
            @test all(isfinite, tendency_x) && all(isfinite, tendency_y)
            @test tendency_x ≈ reference_x rtol = 1.0e-12
            @test tendency_y ≈ reference_y rtol = 1.0e-12

            # a subrange only changes its columns
            tendency_x, tendency_y = diffuse_columns(3:5, couple)
            @test tendency_x[3:5, :] ≈ reference_x[3:5, :] rtol = 1.0e-12
            @test tendency_y[3:5, :] ≈ reference_y[3:5, :] rtol = 1.0e-12
            @test tendency_x[[1, 2, 6, 7], :] == source_x[[1, 2, 6, 7], :]
            @test tendency_y[[1, 2, 6, 7], :] == source_y[[1, 2, 6, 7], :]
        end
    end
end

# boundary layer height of a single column (vectors from top to surface) as defined: where the bulk Richardson
# number relative to the lowermost layer first reaches Ri_c from below, linearly interpolated
function reference_boundary_layer_height(u, v, Tᵥ, Φ, Φₛ, g, cₚ, Ri_c)
    nlayers = length(Tᵥ)
    z(k) = (Φ[k] - Φₛ) / g
    Ri(k) = g * z(k) * (cₚ * Tᵥ[k] + Φ[k] - cₚ * Tᵥ[nlayers] - Φ[nlayers]) / (cₚ * Tᵥ[nlayers] * max(u[k]^2 + v[k]^2, 1))
    for k in (nlayers - 1):-1:1
        if Ri(k) >= Ri_c
            weight = (Ri_c - Ri(k + 1)) / (Ri(k) - Ri(k + 1))
            return z(k + 1) + weight * (z(k) - z(k + 1))
        end
    end
    return z(1)
end

@testset "Boundary layer height" begin
    ncolumns, nlayers = 200, 12
    g, cₚ = 9.81, 1004.0
    rng = Xoshiro(4)
    orography = 2000 * rand(rng, ncolumns)
    z = repeat(reshape(collect(range(15000.0, 100.0, length = nlayers)), 1, :), ncolumns) .* (1 .+ 0.2 * rand(rng, ncolumns))
    Φ = g * (z .+ orography)
    lapse_rate = 4.0e-3 .+ 8.0e-3 * rand(rng, ncolumns)       # [K/m] from stable to unstable (dry adiabatic 9.8e-3)
    Tᵥ = 290 .- lapse_rate .* z .+ 0.3 * randn(rng, ncolumns, nlayers)
    u = 15 * randn(rng, ncolumns, nlayers)
    v = 15 * randn(rng, ncolumns, nlayers)

    for Ri_c in (0.5, 10.0, 1.0e6)
        reference = [reference_boundary_layer_height(u[ij, :], v[ij, :], Tᵥ[ij, :], Φ[ij, :], g * orography[ij], g, cₚ, Ri_c) for ij in 1:ncolumns]
        height = fill(NaN, ncolumns)
        SpeedyWeather.boundary_layer_height!(1:ncolumns, height, u, v, Tᵥ, Φ, orography, g, cₚ, Ri_c)
        @test height ≈ reference rtol = 1.0e-8
        @test all(z[:, end] * (1 - 1.0e-10) .<= height .<= z[:, 1] * (1 + 1.0e-10))     # between lowermost and uppermost layer

        # columns of a subrange only
        height = fill(-1.0, ncolumns)
        SpeedyWeather.boundary_layer_height!(11:20, height, u, v, Tᵥ, Φ, orography, g, cₚ, Ri_c)
        @test height[11:20] ≈ reference[11:20] rtol = 1.0e-8
        @test all(height[[1:10; 21:ncolumns]] .== -1)
    end

    # both cases occur: a crossing inside the column, and none (height of the uppermost layer)
    reference = [reference_boundary_layer_height(u[ij, :], v[ij, :], Tᵥ[ij, :], Φ[ij, :], g * orography[ij], g, cₚ, 10.0) for ij in 1:ncolumns]
    @test 0 < count(reference .< z[:, 1]) < ncolumns || count(reference .< z[:, 1]) == ncolumns
    height = fill(NaN, ncolumns)
    SpeedyWeather.boundary_layer_height!(1:ncolumns, height, u, v, Tᵥ, Φ, orography, g, cₚ, 1.0e6)
    @test height ≈ z[:, 1]
end

# model with vertical diffusion after a day, so that low-level winds and a boundary layer develop
function vertical_diffusion_simulation(Model = PrimitiveWetModel; nlayers = 8, kwargs...)
    spectral_grid = SpectralGrid(truncation = 32; nlayers)
    model = Model(spectral_grid; vertical_diffusion = BulkRichardsonDiffusion(spectral_grid; kwargs...))
    simulation = initialize!(model)
    run!(simulation, period = Day(1))
    return simulation
end

# air density over the distance between layers ρ/Δz [kg/m⁴] at the interfaces k+½ below every layer, from
# the state of the simulation as in vertical_conductance!, zero for the surface interface
function density_over_distance(simulation)
    (; model, variables) = simulation
    (; atmosphere, time_stepping) = model
    (; nlayers, vertical_coordinates) = model.geometry
    diffusion = model.vertical_diffusion
    T = SpeedyWeather.get_prognostic_step(variables.grid.temperature, time_stepping, diffusion).data
    q = haskey(variables.grid, :humidity) ?
        SpeedyWeather.get_prognostic_step(variables.grid.humidity, time_stepping, diffusion).data : zero(T)
    virtual_temperature = SpeedyWeather.virtual_temperature.(T, q, Ref(atmosphere))
    Φ = variables.dynamics.geopotential.data
    pₛ = variables.parameterizations.surface_pressure.data

    ρ_Δz = zeros(eltype(T), size(T))
    for k in 1:(nlayers - 1), ij in axes(T, 1)
        p = SpeedyWeather.pressure_half(k + 1, pₛ[ij], vertical_coordinates)
        Δz = (Φ[ij, k] - Φ[ij, k + 1]) / model.planet.gravity
        Tᵥ = (virtual_temperature[ij, k] + virtual_temperature[ij, k + 1]) / 2
        ρ_Δz[ij, k] = p / (atmosphere.R_dry * Tᵥ * Δz)
    end
    return ρ_Δz
end

# recompute the conductances ρK/Δz on all grid points and return the diffusivities K [m²/s] of momentum
# and heat they correspond to, the boundary layer height and ρ/Δz
function vertical_diffusivity!(simulation)
    (; model, variables) = simulation
    SpeedyWeather.vertical_conductance!(1:model.geometry.npoints, variables, model.vertical_diffusion, model)
    p = variables.parameterizations
    ρ_Δz = density_over_distance(simulation)
    nlayers = size(ρ_Δz, 2)

    # the conductance of the surface interface is zero, as is the diffusivity returned for it
    diffusivity(conductance) =
        [k < nlayers ? conductance[ij, k] / ρ_Δz[ij, k] : conductance[ij, k] for ij in axes(ρ_Δz, 1), k in 1:nlayers]
    K_momentum = diffusivity(p.vertical_conductance_momentum.data)
    K_heat = diffusivity(p.vertical_conductance_heat.data)
    return K_momentum, K_heat, copy(p.boundary_layer_height.data), ρ_Δz
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
            @test all(K .<= max_diffusivity * (1 + 1.0e-5))
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
        K_cold, _, h_cold = diffusivity(-10, 10)        # stable surface layer, 0 < Ri < Ri_cutoff
        K_very_cold, _, _ = diffusivity(-20, 0.5)       # Ri ≥ Ri_cutoff at the surface, no boundary layer profile

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
        _, K_heat, _, ρ_Δz = vertical_diffusivity!(simulation)
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

@testset "BulkRichardsonDiffusion for ranges of columns" begin
    simulation = vertical_diffusion_simulation()
    (; model, variables) = simulation
    (; spectral_grid) = model
    diffusion = model.vertical_diffusion
    npoints = model.geometry.npoints
    p = variables.parameterizations
    names = (:u, :v, :temperature, :humidity)
    tendencies = map(name -> SpeedyWeather.get_tendency_step(getfield(variables.tendencies.grid, name), model.time_stepping, diffusion).data, names)
    initial = map(copy, tendencies)     # of the parameterizations of the last time step, here the source of the diffusion

    # tendencies, conductances and boundary layer height after `run!`, always from the same tendencies
    function diffuse!(run!)
        foreach((tendency, tendency_initial) -> tendency .= tendency_initial, tendencies, initial)
        run!()
        return (;
            tendencies = map(copy, tendencies), conductance_momentum = copy(p.vertical_conductance_momentum.data),
            conductance_heat = copy(p.vertical_conductance_heat.data), height = copy(p.boundary_layer_height.data),
        )
    end
    diffuse!(diffusion::BulkRichardsonDiffusion, columns = 1:npoints) = diffuse!(() -> SpeedyWeather.vertical_diffusion!(columns, variables, diffusion, model))

    @testset "All columns at once, column by column, a subrange" begin
        all_columns = diffuse!(() -> SpeedyWeather.column_parameterization_cpu!(variables, diffusion, model))   # as on CPU
        column_by_column = diffuse!(() -> foreach(ij -> SpeedyWeather.parameterization!(ij, variables, diffusion, model), 1:npoints))   # as in a kernel
        for (a, b) in zip(all_columns.tendencies, column_by_column.tendencies)
            @test all(isfinite, a)
            @test a ≈ b rtol = 1.0e-6
        end
        @test all_columns.conductance_momentum ≈ column_by_column.conductance_momentum rtol = 1.0e-6
        @test all_columns.conductance_heat ≈ column_by_column.conductance_heat rtol = 1.0e-6
        @test all_columns.height ≈ column_by_column.height rtol = 1.0e-6
        @test all(map((a, b) -> any(a .!= b), all_columns.tendencies, initial))     # every variable is diffused

        columns = 101:200
        others = setdiff(1:npoints, columns)
        subrange = diffuse!(diffusion, columns)
        for (a, b, c) in zip(subrange.tendencies, all_columns.tendencies, initial)
            @test a[columns, :] ≈ b[columns, :] rtol = 1.0e-6
            @test a[others, :] == c[others, :]
        end
    end

    @testset "Options: coupled = $couple_surface_fluxes" for couple_surface_fluxes in (true, false)
        reference = diffuse!(BulkRichardsonDiffusion(spectral_grid; couple_surface_fluxes))
        for diffuse_momentum in (true, false), diffuse_static_energy in (true, false), diffuse_humidity in (true, false)
            options = BulkRichardsonDiffusion(spectral_grid; couple_surface_fluxes, diffuse_momentum, diffuse_static_energy, diffuse_humidity)
            result = diffuse!(options)
            diffused = (diffuse_momentum, diffuse_momentum, diffuse_static_energy, diffuse_humidity)
            for (i, is_diffused) in enumerate(diffused)
                if is_diffused      # as with all variables diffused, also if solved alone
                    @test result.tendencies[i] ≈ reference.tendencies[i] rtol = 1.0e-6
                else                # untouched
                    @test result.tendencies[i] == initial[i]
                end
            end
        end
    end

    @testset "Surface Richardson number stored by the drag" begin
        drag = model.boundary_layer.drag
        @test drag isa BulkRichardsonDrag
        g = model.planet.gravity
        Φ = variables.dynamics.geopotential.data
        orography = model.orography.orography.data
        foreach(ij -> SpeedyWeather.boundary_layer_drag!(ij, variables, drag, model), 1:npoints)

        # as the diffusion would compute it without the stored one
        recomputed = map(1:npoints) do ij
            ΔΦ₀ = Φ[ij, end] - g * orography[ij]
            Ri_land, Ri_ocean = SpeedyWeather.bulk_richardson_surface(ij, ΔΦ₀, variables, model.atmosphere, diffusion, model.time_stepping)
            land_fraction = model.land_sea_mask.land_fraction[ij]
            land_fraction * Ri_land + (1 - land_fraction) * Ri_ocean
        end
        stored = p.surface_bulk_richardson.data
        @test all(isfinite, stored)
        @test all(isapprox.(stored, recomputed, rtol = 1.0e-3, atol = 1.0e-3))
        @test SpeedyWeather.surface_richardson(17, 0.0f0, variables, diffusion, model) === stored[17]
    end

    @testset "Drag without a stored surface Richardson number" begin
        boundary_layer = BoundaryLayer(spectral_grid; drag = ConstantDrag(spectral_grid))
        model_constant_drag = PrimitiveWetModel(spectral_grid; boundary_layer)
        simulation_constant_drag = initialize!(model_constant_drag)
        run!(simulation_constant_drag, period = Day(1))
        parameterizations = simulation_constant_drag.variables.parameterizations
        @test !haskey(parameterizations, :surface_bulk_richardson)
        @test !simulation_constant_drag.model.feedback.nans_detected
        @test all(isfinite, parameterizations.vertical_conductance_momentum)
        @test maximum(parameterizations.vertical_conductance_heat) > 0

        # then computed by the diffusion itself
        ΔΦ₀ = 500.0f0
        Ri = SpeedyWeather.surface_richardson(1, ΔΦ₀, simulation_constant_drag.variables, model_constant_drag.vertical_diffusion, model_constant_drag)
        @test isfinite(Ri)
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
        (; vertical_conductance_momentum, vertical_conductance_heat, boundary_layer_height) = simulation.variables.parameterizations
        @test all(isfinite, vertical_conductance_momentum)
        @test all(isfinite, vertical_conductance_heat)
        @test maximum(vertical_conductance_heat) > 0
        @test all(boundary_layer_height .> 0)
    end

    # disabled diffusion leaves the tendencies unchanged
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    diffusion = BulkRichardsonDiffusion(spectral_grid; diffuse_momentum = false, diffuse_static_energy = false, diffuse_humidity = false)
    model = PrimitiveWetModel(spectral_grid; vertical_diffusion = diffusion)
    simulation = initialize!(model)
    run!(simulation, steps = 3)
    @test all(simulation.variables.parameterizations.vertical_conductance_heat .== 0)   # never computed

    # friction velocity and surface Richardson number from the boundary layer, both zero without its
    # drag and surface wind, which switches the boundary-layer profile of the diffusivity off
    (; variables) = simulation
    ΔΦ₀ = 500.0f0       # geopotential of the lowermost layer relative to the surface
    drag, surface_wind_speed = SpeedyWeather.boundary_layer_drag_and_wind(variables.parameterizations)
    @test SpeedyWeather.friction_velocity(drag, surface_wind_speed, 1) > 0
    @test isfinite(SpeedyWeather.surface_richardson(1, ΔΦ₀, variables, diffusion, model))
    no_boundary_layer = (; parameterizations = (; variables.parameterizations.vertical_conductance_momentum))
    @test SpeedyWeather.boundary_layer_drag_and_wind(no_boundary_layer.parameterizations) === (nothing, nothing)
    @test SpeedyWeather.friction_velocity(nothing, nothing, 1) == 0
    @test SpeedyWeather.surface_richardson(1, ΔΦ₀, no_boundary_layer, diffusion, model) === 0.0f0
end
