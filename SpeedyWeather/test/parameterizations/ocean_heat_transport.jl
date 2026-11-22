using Random

# a simulation with a slab ocean and diffusion, at truncation, grid and number format
function heat_diffusion_simulation(; truncation = 32, Grid = FullGaussianGrid, NF = Float64, kwargs...)
    spectral_grid = SpectralGrid(; truncation, Grid, NF, nlayers = 4)
    ocean = SlabOcean(spectral_grid, heat_transport = OceanHeatDiffusion(spectral_grid; kwargs...))
    model = PrimitiveWetModel(spectral_grid; ocean)
    return initialize!(model)
end

# SST tendency from the heat transport only, for sea surface temperature `T`
function heat_transport_tendency(simulation, T)
    (; model, variables) = simulation
    sst = SpeedyWeather.get_prognostic_step(variables.prognostic.ocean.sea_surface_temperature, model.time_stepping.ocean, model.ocean)
    dsst = SpeedyWeather.get_tendency_step(variables.tendencies.ocean.sea_surface_temperature, model.time_stepping.ocean, model.ocean)
    sst .= T
    dsst .= 0
    SpeedyWeather.heat_transport!(dsst, variables, model.ocean.heat_transport, model.ocean, model)
    return copy(dsst)
end

@testset "OceanHeatDiffusion" begin
    Random.seed!(1)

    @testset "Default SlabOcean has no heat transport" begin
        sg = SpectralGrid(truncation = 21, nlayers = 4)
        @test SlabOcean(sg).heat_transport isa NoOceanHeatTransport
        model = PrimitiveWetModel(sg, ocean = SlabOcean(sg))
        simulation = initialize!(model)
        run!(simulation, steps = 2)
        sst_default = copy(simulation.variables.prognostic.ocean.sea_surface_temperature)

        # same with an explicit no-op, and the default run is reproducible
        model = PrimitiveWetModel(sg, ocean = SlabOcean(sg, heat_transport = NoOceanHeatTransport()))
        simulation = initialize!(model)
        run!(simulation, steps = 2)
        @test sst_default == simulation.variables.prognostic.ocean.sea_surface_temperature
        @test !any(isnan, sst_default)
    end

    for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
        @testset "Conservation, uniformity and coasts, $Grid" begin
            simulation = heat_diffusion_simulation(; Grid)
            (; model) = simulation
            (; land_fraction) = model.land_sea_mask
            ht = model.ocean.heat_transport
            land = land_fraction.data .>= 1
            ocean_fraction = 1 .- land_fraction.data
            @test any(land) && any(.!land)          # real land-sea mask
            @test ht.diffusivity_used ≈ 6250

            # rows sum to zero, heat-weighted columns sum to zero
            @test all(abs.(sum(ht.weights, dims = 2)) .< 1.0e-12 * ht.max_rate)

            # uniform field stays uniform, with land
            @test maximum(abs, heat_transport_tendency(simulation, 285.0)) < 1.0e-12 * ht.max_rate * 285

            # heat conserved to round-off with random temperatures
            for _ in 1:3
                T = 280 .+ 10 .* rand(length(land_fraction.data))
                dT = heat_transport_tendency(simulation, T)
                heat = sum(ht.area .* ocean_fraction .* dT.data)
                scale = sum(ht.area .* ocean_fraction .* abs.(dT.data))
                @test abs(heat) < 1.0e-12 * scale
                @test any(dT.data .!= 0)
            end

            # no flux into land: land points are untouched, and the ocean does not see land temperatures
            T = 280 .+ 10 .* rand(length(land_fraction.data))
            dT = heat_transport_tendency(simulation, T)
            @test all(dT.data[land] .== 0)
            T2 = copy(T)
            T2[land] .+= 100
            @test heat_transport_tendency(simulation, T2).data ≈ dT.data

            # conservation over many steps of explicit Euler with random initial temperatures
            Δt = model.time_stepping.Δt
            sst = SpeedyWeather.get_prognostic_step(simulation.variables.prognostic.ocean.sea_surface_temperature, model.time_stepping.ocean, model.ocean)
            weights = ht.area .* ocean_fraction
            T = 280 .+ 10 .* rand(length(land_fraction.data))
            heat0 = sum(weights[.!land] .* T[.!land])
            for _ in 1:500
                dT = heat_transport_tendency(simulation, T)
                T = T .+ Δt .* dT.data
            end
            @test abs(sum(weights[.!land] .* T[.!land]) - heat0) < 1.0e-12 * abs(heat0)
            @test extrema(T[.!land])[2] - extrema(T[.!land])[1] < 10    # diffused
        end
    end

    @testset "Laplacian accuracy" begin
        # spherical harmonic l=3, m=2 (away from the poles): ∇²Y = -12 Y / R²
        Y(λ, φ) = cos(φ)^2 * sin(φ) * cos(2λ)
        for (truncation, Grid, tolerance) in ((32, FullGaussianGrid, 0.02), (32, OctahedralGaussianGrid, 0.02), (64, FullGaussianGrid, 0.01))
            simulation = heat_diffusion_simulation(; truncation, Grid, diffusivity = 1.0)
            (; model) = simulation
            model.land_sea_mask.land_fraction .= 0                      # aquaplanet
            initialize!(model.ocean.heat_transport, model.ocean, model)
            londs, latds = RingGrids.get_londlatds(model.spectral_grid.grid)
            T = Y.(deg2rad.(londs), deg2rad.(latds))
            dT = heat_transport_tendency(simulation, T)
            exact = -12 / model.geometry.radius[]^2 .* T                # D = 1 m²/s
            area = model.ocean.heat_transport.area      # area-weighted: the polar caps are tiny but have the largest errors
            relative_error = sqrt(sum(area .* abs2.(dT.data .- exact)) / sum(area .* abs2.(exact)))
            @info "Laplacian relative error" Grid truncation relative_error
            @test relative_error < tolerance
        end
    end

    @testset "Stability number" begin
        for truncation in (32, 42)
            simulation = heat_diffusion_simulation(; truncation)
            (; model) = simulation
            number = SpeedyWeather.stability_number(model.ocean.heat_transport, model.time_stepping.Δt)
            @info "Stability number" truncation number
            @test 0 < number < 1
        end
        @test_logs (:warn, r"unstable") match_mode = :any heat_diffusion_simulation(; diffusivity = 1.0e7)
    end

    @testset "Runs with sea ice" begin
        simulation = heat_diffusion_simulation()
        run!(simulation, steps = 10)
        @test simulation.model.feedback.nans_detected == false
        @test !any(isnan, simulation.variables.prognostic.ocean.sea_surface_temperature)
    end
end
