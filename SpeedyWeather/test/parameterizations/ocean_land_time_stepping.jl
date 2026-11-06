@testset "EulerForward" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    euler = EulerForward(spectral_grid)
    model = PrimitiveWetModel(spectral_grid)
    clock = Clock()

    var = zeros(spectral_grid.grid, 2)
    var .= 280
    tendency = zeros(spectral_grid.grid, 2)
    tendency .= 1.0e-4
    SpeedyWeather.update_prognostic!(var, tendency, clock, euler, nothing, model)
    @test all(var .≈ 280 + euler.Δt * 1.0e-4)

    var .= 280      # scale divides the time step
    SpeedyWeather.update_prognostic!(var, tendency, clock, euler, nothing, model, 2)
    @test all(var .≈ 280 + euler.Δt / 2 * 1.0e-4)

    @test SpeedyWeather.prognostic_grid_steps(euler, model) == 1
    @test SpeedyWeather.tendency_grid_steps(euler, model) == 1
end

@testset "Ocean and land time stepping with leapfrog" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    (; variables) = simulation
    (; clock) = variables.prognostic
    leapfrog = model.time_stepping

    # default EulerForward for ocean and land, with the leapfrog's Δt
    @test SpeedyWeather.time_stepper(model.time_stepping, :ocean) === leapfrog.ocean
    @test SpeedyWeather.time_stepper(model.time_stepping, :land) === leapfrog.land
    @test leapfrog.ocean isa EulerForward
    @test leapfrog.land isa EulerForward
    @test leapfrog.land.Δt == leapfrog.ocean.Δt == leapfrog.Δt

    # set! also changes the time step of ocean and land
    set!(model, Δt = Minute(10))
    @test leapfrog.land.Δt == leapfrog.ocean.Δt == leapfrog.Δt == 600

    # Δt/2 on the first Euler + first leapfrog step as the clock only advances by Δt/2 each, Δt thereafter
    Δt = leapfrog.Δt
    for (step_counter, expected) in ((0, Δt / 2), (1, Δt / 2), (2, Δt), (100, Δt))
        clock.step_counter = step_counter
        @test SpeedyWeather.time_step(model.time_stepping, :land, clock) == expected
        @test SpeedyWeather.time_step(model.time_stepping, :ocean, clock) == expected
    end

    # only 1 step allocated for ocean and land prognostic variables and tendencies
    for var in (
            variables.prognostic.land.soil_temperature, variables.tendencies.land.soil_temperature,
            variables.prognostic.ocean.sea_surface_temperature, variables.tendencies.ocean.sea_surface_temperature,
        )
        @test SpeedyWeather.nsteps(var) == 1
    end
    run!(simulation, steps = 4)
    @test simulation.model.feedback.nans_detected == false

    # NCycleLorenz uses itself for ocean and land by default
    ncycle = NCycleLorenz(spectral_grid)
    @test SpeedyWeather.time_stepper(ncycle, :land) === ncycle
    @test SpeedyWeather.time_stepper(ncycle, :ocean) === ncycle
end

@testset "NCycleLorenz for ocean and land with leapfrog" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    time_stepping = Leapfrog(spectral_grid, ocean = NCycleLorenz(spectral_grid), land = NCycleLorenz(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    (; variables) = simulation
    @test time_stepping.land.Δt == time_stepping.Δt
    # the N-cycle decides on the steps: 1 prognostic step, 2 tendency steps (F, G)
    @test SpeedyWeather.nsteps(variables.prognostic.land.soil_temperature) == 1
    @test SpeedyWeather.nsteps(variables.tendencies.land.soil_temperature) == 2
    @test SpeedyWeather.nsteps(variables.tendencies.ocean.sea_surface_temperature) == 2
    run!(simulation, steps = 4)
    @test simulation.model.feedback.nans_detected == false
end

@testset "EulerForward for ocean and land with NCycleLorenz" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    time_stepping = NCycleLorenz(spectral_grid, ocean = EulerForward(spectral_grid), land = EulerForward(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    @test SpeedyWeather.time_stepper(model.time_stepping, :land) === time_stepping.land
    @test time_stepping.land.Δt == time_stepping.ocean.Δt == time_stepping.Δt
    @test SpeedyWeather.nsteps(simulation.variables.tendencies.land.soil_temperature) == 1
    # no run! as NCycleLorenz does not support primitive equation models yet
end

@testset "Leapfrogged ocean and land" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    time_stepping = Leapfrog(spectral_grid, ocean = nothing, land = nothing)
    model = PrimitiveWetModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    @test SpeedyWeather.time_stepper(model.time_stepping, :land) === time_stepping
    @test SpeedyWeather.nsteps(simulation.variables.prognostic.land.soil_temperature) == 2
    @test SpeedyWeather.nsteps(simulation.variables.prognostic.ocean.sea_surface_temperature) == 2
    run!(simulation, steps = 4)
    @test simulation.model.feedback.nans_detected == false
end

@testset "Snow time stepped by the land time stepper" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    for land_time_stepper in (EulerForward(spectral_grid), nothing)    # nothing = leapfrogged
        time_stepping = Leapfrog(spectral_grid, land = land_time_stepper)
        model = PrimitiveWetModel(spectral_grid; time_stepping)
        @test model.land.snow isa SnowModel
        simulation = initialize!(model)
        (; variables) = simulation
        snow_depth = variables.prognostic.land.snow_depth
        nsteps = isnothing(land_time_stepper) ? 2 : 1
        @test SpeedyWeather.nsteps(snow_depth) == nsteps
        @test haskey(variables.tendencies.land, :snow_depth)

        # filter! keeps snow depth in [0, snow_depth_cap] for all steps
        snow_depth .= 1000
        snow_depth[1, :] .= -1
        filter!(variables, model.land.snow, model)
        @test all(snow_depth[1, :] .== 0)
        @test all(snow_depth[2:end, :] .== model.land.snow.snow_depth_cap)

        snow_depth .= 0
        run!(simulation, steps = 4)
        @test simulation.model.feedback.nans_detected == false
        @test all(0 .<= snow_depth .<= model.land.snow.snow_depth_cap)
    end
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

# record after every time step whether ocean points are below freezing and the total sea ice
mutable struct SeaIceLog <: SpeedyWeather.AbstractCallback
    ncold::Vector{Int}
    ice::Vector{Float64}
end
SpeedyWeather.initialize!(::SeaIceLog, args...) = nothing
SpeedyWeather.finalize!(::SeaIceLog, args...) = nothing
function SpeedyWeather.callback!(log::SeaIceLog, vars, model)
    sst = vars.prognostic.ocean.sea_surface_temperature
    ocean = model.land_sea_mask.land_fraction.data .< 1
    push!(log.ncold, count(ocean .& (sst.data .< model.sea_ice.freezing_temperature - 1.0e-3)))
    push!(log.ice, sum(vars.prognostic.ocean.sea_ice_concentration))
    return nothing
end

@testset "Sea ice freezing and snow melt cap independent of the time stepper" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    for time_stepper in (EulerForward(spectral_grid), NCycleLorenz(spectral_grid))
        time_stepping = Leapfrog(spectral_grid, ocean = time_stepper, land = time_stepper)
        model = PrimitiveWetModel(spectral_grid; time_stepping)
        log = SeaIceLog(Int[], Float64[])
        add!(model, log)
        simulation = initialize!(model)
        run!(simulation, steps = 8)
        @test simulation.model.feedback.nans_detected == false
        @test all(log.ncold .== 0)          # SST restored to freezing after every step
        @test issorted(log.ice)             # sea ice grows monotonically in the first steps (no N-cycle flip)
        @test all(0 .<= simulation.variables.prognostic.land.snow_depth .<= model.land.snow.snow_depth_cap)
    end
end
