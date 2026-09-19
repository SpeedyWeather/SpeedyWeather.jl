@testset "Surface Euler time step with leapfrog" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    (; variables) = simulation
    (; clock) = variables.prognostic
    leapfrog = model.time_stepping
    Δt = leapfrog.Δt

    # Δt/2 on the first Euler + first leapfrog step as the clock only advances by Δt/2 each, Δt thereafter
    for (step_counter, expected) in ((0, Δt / 2), (1, Δt / 2), (2, Δt), (100, Δt))
        clock.step_counter = step_counter
        @test SpeedyWeather.surface_time_step(leapfrog, clock) == expected
    end

    # Euler forward from the 1st step, written into both steps
    soil_temperature = variables.prognostic.land.soil_temperature
    soil_temperature_tendency = variables.tendencies.land.soil_temperature
    @test ArrayDimensions.hastime(soil_temperature)

    get_step(soil_temperature, 1) .= 280
    get_step(soil_temperature, 2) .= 300     # should be overwritten
    soil_temperature_tendency .= 1.0e-4

    clock.step_counter = 2
    SpeedyWeather.update_prognostic_surface!(soil_temperature, soil_temperature_tendency, clock, leapfrog, model.implicit, model)
    @test all(get_step(soil_temperature, 1) .≈ 280 + Δt * 1.0e-4)
    @test get_step(soil_temperature, 1) == get_step(soil_temperature, 2)

    # other time steppers fall back to their update_prognostic!
    @test SpeedyWeather.surface_time_step(NCycleLorenz(spectral_grid), clock) == NCycleLorenz(spectral_grid).Δt
end

@testset "LandBucketMoisture excess water is conserved" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    (; variables) = simulation
    soil = model.land.soil_moisture
    @test soil isa LandBucketMoisture

    land_fraction = Array(model.land_sea_mask.land_fraction.data)
    ij = findfirst(>(0), land_fraction)

    γ = model.land.thermodynamics.field_capacity
    f₁ = γ * model.land.geometry.layer_thickness[1]
    f₂ = γ * model.land.geometry.layer_thickness[2]
    p = soil.infiltration_fraction

    soil_moisture = variables.prognostic.land.soil_moisture
    R = variables.parameterizations.land.river_runoff
    soil_moisture[ij, 1, :] .= 1.2     # 0.2 excess in top layer
    soil_moisture[ij, 2, :] .= 0.5
    R[ij] = 0

    water_before = f₁ * 1.2 + f₂ * 0.5
    filter!(variables, soil, model)

    for step in 1:size(soil_moisture, 3)
        @test soil_moisture[ij, 1, step] ≈ 1
        @test soil_moisture[ij, 2, step] ≈ 0.5 + p * 0.2 * f₁ / f₂
    end
    @test R[ij] ≈ (1 - p) * 0.2 * f₁                   # runoff only accumulated once
    @test f₁ * soil_moisture[ij, 1, 1] + f₂ * soil_moisture[ij, 2, 1] + R[ij] ≈ water_before

    # everything in [0, 1]
    @test all(0 .<= soil_moisture .<= 1)
end

@testset "Ocean, sea ice, land type hierarchy" begin
    @test LandBucketTemperature <: SpeedyWeather.AbstractLandTemperature
    @test SeasonalLandTemperature <: SpeedyWeather.AbstractLandTemperature
    @test LandBucketMoisture <: SpeedyWeather.AbstractSoilMoisture
    @test SeasonalSoilMoisture <: SpeedyWeather.AbstractSoilMoisture
    @test SnowModel <: SpeedyWeather.AbstractSnow
    @test ThermodynamicSeaIce <: SpeedyWeather.AbstractSeaIce
    @test PrescribedSeaIce <: SpeedyWeather.AbstractSeaIce
    @test SlabOcean <: SpeedyWeather.AbstractOcean
    @test all(T -> T <: SpeedyWeather.AbstractLandComponent, (LandBucketTemperature, LandBucketMoisture, SnowModel))
end

# a custom ocean that only subtypes AbstractOcean still gets its sea surface temperature
struct TestCustomOcean <: SpeedyWeather.AbstractOcean end
SpeedyWeather.initialize!(::TestCustomOcean, ::PrimitiveEquation) = nothing
SpeedyWeather.initialize!(vars::Variables, ::TestCustomOcean, ::PrimitiveEquation) =
    (vars.prognostic.ocean.sea_surface_temperature .= 290; nothing)
SpeedyWeather.timestep!(::Variables, ::TestCustomOcean, ::PrimitiveEquation) = nothing

@testset "Custom ocean and prescribed land variables" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid; ocean = TestCustomOcean())
    simulation = initialize!(model)
    @test haskey(simulation.variables.prognostic.ocean, :sea_surface_temperature)
    @test !ArrayDimensions.hastime(simulation.variables.prognostic.ocean.sea_surface_temperature)

    # snow does not define soil temperature, so a prescribed land temperature stays without time dimension
    land = LandModel(spectral_grid; temperature = SeasonalLandTemperature(spectral_grid), snow = SnowModel(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; land)
    simulation = initialize!(model)
    @test !ArrayDimensions.hastime(simulation.variables.prognostic.land.soil_temperature)
    run!(simulation, steps = 4)
    @test simulation.model.feedback.nans_detected == false
end

@testset "Latent heat of fusion for snow melt" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    @test SpeedyWeather.latent_heat_fusion(PrimitiveWetModel(spectral_grid).atmosphere) ≈ 3.3e5
    @test SpeedyWeather.latent_heat_fusion(PrimitiveDryModel(spectral_grid).atmosphere) == 0
end
