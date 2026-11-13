@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))

        initialize!(model.radiation, model)

        vars = Variables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.radiation.longwave, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    @testset for T in (FriersonLongwaveTransmissivity, ByrneOGormanLongwaveTransmissivity, TransparentLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
        @test model.radiation.longwave === longwave
        initialize!(model.radiation, model)

        vars = Variables(model)

        # transmissivity depends on pressure thickness and thereofre surface pressure should be nonzero
        vars.parameterizations.surface_pressure .= 1.0e5
        t = SpeedyWeather.transmissivity!(1, vars, model.radiation.longwave.transmissivity, model)
        for ij in 2:model.spectral_grid.npoints
            SpeedyWeather.transmissivity!(ij, vars, model.radiation.longwave.transmissivity, model)
        end

        @test all(0 .< t .<= 1)
    end
end

@testset "Byrne and O'Gorman humidity-dependent longwave transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    transmissivity = ByrneOGormanLongwaveTransmissivity(spectral_grid)
    (; dry_absorption, water_vapor_absorption, co2_sensitivity, reference_co2) = transmissivity
    longwave = OneBandLongwave(spectral_grid; transmissivity)
    model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
    vars = initialize!(model).variables     # initialize! sets CO₂ from model.greenhouse_gases
    vars.parameterizations.surface_pressure .= 1.0e5

    # dry: only the dry optical depth over the column, at the default CO₂
    co2 = vars.prognostic.greenhouse_gases.co2[]
    @test co2 == SpeedyWeather.DEFAULT_CO2
    τ_dry = dry_absorption + co2_sensitivity * log(co2 / reference_co2)
    vars.grid.humidity .= 0
    t_dry = copy(SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :])
    @test prod(t_dry) ≈ exp(-τ_dry) rtol = 1.0e-5

    # moist: lower transmissivity in every layer
    vars.grid.humidity .= 0.005
    t_moist = copy(SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :])
    @test all(t_moist .< t_dry)
    @test prod(t_moist) ≈ exp(-(τ_dry + water_vapor_absorption * 0.005)) rtol = 1.0e-4

    # doubling CO₂ adds co2_sensitivity * log(2) to the column optical depth
    vars.prognostic.greenhouse_gases.co2[] = 2co2
    t_2xco2 = SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :]
    @test all(t_2xco2 .< t_moist)
    @test prod(t_2xco2) ≈ prod(t_moist) * exp(-co2_sensitivity * log(2)) rtol = 1.0e-4

    # default for the wet model
    @test PrimitiveWetModel(spectral_grid).radiation.longwave.transmissivity isa ByrneOGormanLongwaveTransmissivity

    # without greenhouse gases: CO₂ at the reference concentration
    model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave), greenhouse_gases = nothing)
    vars = initialize!(model).variables
    vars.parameterizations.surface_pressure .= 1.0e5
    vars.grid.humidity .= 0
    t_ref = SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :]
    @test prod(t_ref) ≈ exp(-dry_absorption) rtol = 1.0e-5
end 

@testset "Stratospheric longwave emission" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    initialize!(model.radiation.longwave, model)
    @test model.radiation.longwave.radiative_transfer.stratospheric_emissivity > 0     # on for wet model
    @test PrimitiveDryModel(spectral_grid).radiation.longwave.radiative_transfer.stratospheric_emissivity == 0

    # compare with and without stratospheric emission for the same state
    function column(ϵ)
        radiative_transfer = OneBandLongwaveRadiativeTransfer(spectral_grid, stratospheric_emissivity = ϵ)
        longwave = OneBandLongwave(spectral_grid; radiative_transfer)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
        simulation = initialize!(model)
        vars = simulation.variables
        vars.parameterizations.surface_pressure .= 1.0e5
        vars.grid.temperature .= 250
        vars.tendencies.grid.temperature .= 0
        ij = 1
        SpeedyWeather.parameterization!(ij, vars, model.radiation.longwave, model)
        dTdt = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, model.radiation.longwave)
        return Array(dTdt[ij, :]), vars.parameterizations.outgoing_longwave[ij], model
    end

    dTdt0, olr0, model = column(0)
    dTdt1, olr1, _ = column(0.05)
    coordinates = model.geometry.vertical_coordinates
    stratosphere = [SpeedyWeather.pressure_above(k, 1.0e5, coordinates) < 14000 for k in 1:spectral_grid.nlayers]

    @test all(dTdt1[stratosphere] .< dTdt0[stratosphere])          # cools the stratosphere
    @test dTdt1[.!stratosphere] ≈ dTdt0[.!stratosphere]             # troposphere unchanged
    @test olr1 > olr0                                               # emitted to space

    # energy conservation: extra OLR equals column-integrated extra cooling
    cₚ = model.atmosphere.heat_capacity
    g = model.planet.gravity
    Δp = [SpeedyWeather.pressure_thickness(k, 1.0e5, coordinates) for k in 1:spectral_grid.nlayers]
    @test sum((dTdt0 .- dTdt1) .* Δp) * cₚ / g ≈ olr1 - olr0 rtol = 1.0e-3
end
