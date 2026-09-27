@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave_radiation = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)

        initialize!(model.longwave_radiation, model)

        vars = Variables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    @testset for T in (FriersonLongwaveTransmissivity, TransparentLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave_radiation = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)
        initialize!(model.longwave_radiation, model)

        vars = Variables(model)

        # transmissivity depends on pressure thickness and thereofre surface pressure should be nonzero
        vars.parameterizations.surface_pressure .= 1.0e5
        t = SpeedyWeather.transmissivity!(1, vars, model.longwave_radiation.transmissivity, model)
        for ij in 2:model.spectral_grid.npoints
            SpeedyWeather.transmissivity!(ij, vars, model.longwave_radiation.transmissivity, model)
        end

        @test all(0 .< t .<= 1)
    end
end

@testset "Stratospheric longwave emission" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    initialize!(model.longwave_radiation, model)
    @test model.longwave_radiation.radiative_transfer.stratospheric_emissivity > 0     # on for wet model
    @test PrimitiveDryModel(spectral_grid).longwave_radiation.radiative_transfer.stratospheric_emissivity == 0

    # compare with and without stratospheric emission for the same state
    function column(ϵ)
        radiative_transfer = OneBandLongwaveRadiativeTransfer(spectral_grid, stratospheric_emissivity = ϵ)
        longwave_radiation = OneBandLongwave(spectral_grid; radiative_transfer)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)
        simulation = initialize!(model)
        vars = simulation.variables
        vars.parameterizations.surface_pressure .= 1.0e5
        vars.grid.temperature .= 250
        vars.tendencies.grid.temperature .= 0
        ij = 1
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
        dTdt = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, model.longwave_radiation)
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
