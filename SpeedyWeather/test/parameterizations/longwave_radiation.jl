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

@testset "Surface longwave emission over sea ice" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for LW in (JeevanjeeRadiation, OneBandLongwave)
        longwave = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
        simulation = initialize!(model)
        vars = simulation.variables
        @test haskey(vars.prognostic.ocean, :sea_ice_concentration)

        ϵ = model.radiation.longwave isa JeevanjeeRadiation ? model.radiation.longwave.emissivity_ocean :
            model.radiation.longwave.radiative_transfer.emissivity_ocean
        σ = model.atmosphere.stefan_boltzmann
        ij = 1
        function surface_longwave_up(sea_ice_concentration, surface_air_temperature; sea_surface_temperature = 271)
            vars.parameterizations.surface_pressure .= 1.0e5
            vars.grid.temperature .= 250
            vars.prognostic.ocean.sea_surface_temperature .= sea_surface_temperature
            vars.prognostic.ocean.sea_ice_concentration .= sea_ice_concentration
            vars.parameterizations.surface_air_temperature .= surface_air_temperature
            SpeedyWeather.parameterization!(ij, vars, model.radiation.longwave, model)
            return vars.parameterizations.ocean.surface_longwave_up[ij]
        end

        @test surface_longwave_up(0, 243) ≈ ϵ * σ * 271^4       # open water emits at the sea surface temperature
        @test surface_longwave_up(1, 243) ≈ ϵ * σ * 243^4       # sea ice at the (colder) surface air temperature
        @test surface_longwave_up(1, 280) ≈ ϵ * σ * 271^4       # but never warmer than the water underneath
        @test surface_longwave_up(0.5, 243) ≈ ϵ * σ * (271 / 2 + 243 / 2)^4

        # no ice: exactly the open-water emission, for any surface air temperature
        @test surface_longwave_up(0, 200) == ϵ * σ * 271^4
        @test surface_longwave_up(0, 300) == ϵ * σ * 271^4

        # bounds, monotonicity and number format of the emission temperature
        sst = vars.prognostic.ocean.sea_surface_temperature
        NF = eltype(sst)
        T_eff(ℵ, Ta, S) = (
            sst .= S; vars.prognostic.ocean.sea_ice_concentration .= ℵ;
            vars.parameterizations.surface_air_temperature .= Ta;
            SpeedyWeather.ocean_surface_emission_temperature(ij, vars, sst, model.radiation.longwave, model)
        )
        @test T_eff(0.5, 243, 271) isa NF
        for ℵ in (0, 0.25, 0.5, 1), Ta in (230, 260, 271, 290), S in (271, 280)
            T = T_eff(ℵ, Ta, S)
            @test min(S, Ta) - 1.0e-3 <= T <= S + 1.0e-3                     # between the ice and the water
            @test T_eff(min(ℵ + 0.25, 1), Ta, S) <= T + 1.0e-3                # more ice: not warmer
            @test T_eff(ℵ, Ta + 5, S) >= T - 1.0e-3                           # warmer air: not colder
            @test T_eff(ℵ, Ta + 5, S) - T <= ℵ * 5 + 1.0e-3                   # ... by at most ℵ ΔTa
        end

        # land isolation and energy closure: the land emission does not depend on the sea ice, and the upward
        # longwave the atmosphere receives is the land-sea weighted sum of what land and ocean models lose
        land_fraction = model.land_sea_mask.land_fraction
        function all_columns(ℵ)
            vars.parameterizations.surface_pressure .= 1.0e5
            vars.grid.temperature .= 250
            sst .= 271
            vars.prognostic.land.soil_temperature .= 280
            vars.prognostic.ocean.sea_ice_concentration .= ℵ
            vars.parameterizations.surface_air_temperature .= 243
            for ij in eachindex(land_fraction)
                SpeedyWeather.parameterization!(ij, vars, model.radiation.longwave, model)
            end
            return copy(vars.parameterizations.land.surface_longwave_up.data),
                copy(vars.parameterizations.ocean.surface_longwave_up.data),
                copy(vars.parameterizations.surface_longwave_up.data)
        end
        land_open, ocean_open, total_open = all_columns(0)
        land_ice, ocean_ice, total_ice = all_columns(1)
        @test land_ice == land_open
        @test all(ocean_ice .< ocean_open)
        lf = land_fraction.data
        @test total_ice ≈ lf .* land_ice .+ (1 .- lf) .* ocean_ice
        @test total_open ≈ lf .* land_open .+ (1 .- lf) .* ocean_open

        # Jensen gap (known limitation): the linear blend underestimates the cell-average flux
        # εσ[(1-ℵ)SST⁴ + ℵTa⁴] by ≈ 6εσT̄²ℵ(1-ℵ)ΔT² (second order, T̄ the mean of SST and Ta)
        for (ℵ, Ta) in ((0.5, 241), (0.5, 211), (0.3, 261), (0.8, 231))
            gap = ϵ * σ * ((1 - ℵ) * 271^4 + ℵ * Ta^4) - surface_longwave_up(ℵ, Ta)
            @test gap >= 0
            @test gap ≈ 6 * ϵ * σ * ((271 + Ta) / 2)^2 * ℵ * (1 - ℵ) * (271 - Ta)^2 rtol = 0.05
        end
    end

    # the emission temperature dispatches on the sea ice model
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    vars = simulation.variables
    sst = vars.prognostic.ocean.sea_surface_temperature
    sst .= 271
    vars.prognostic.ocean.sea_ice_concentration .= 1
    vars.parameterizations.surface_air_temperature .= 243
    @test model.sea_ice isa ThermodynamicSeaIce
    @test SpeedyWeather.ocean_surface_emission_temperature(1, vars, sst, model.radiation.longwave, model) ≈ 243
    @test SpeedyWeather.ocean_surface_emission_temperature(1, vars, sst, PrescribedSeaIce(), model.radiation.longwave, model) ≈ 243
    # sea_ice = nothing with a manually set concentration also blends
    @test SpeedyWeather.ocean_surface_emission_temperature(1, vars, sst, nothing, model.radiation.longwave, model) ≈ 243

    # PrescribedSeaIce: the concentration has no time step dimension
    model = PrimitiveWetModel(spectral_grid; sea_ice = PrescribedSeaIce())
    simulation = initialize!(model)
    vars = simulation.variables
    sst = vars.prognostic.ocean.sea_surface_temperature
    sst .= 271
    vars.prognostic.ocean.sea_ice_concentration .= 1
    vars.parameterizations.surface_air_temperature .= 243
    @test SpeedyWeather.ocean_surface_emission_temperature(1, vars, sst, model.radiation.longwave, model) ≈ 243

    # the surface air temperature is computed (boundary layer) before the radiation uses it
    @test findfirst(==(:boundary_layer), model.parameterizations) < findfirst(==(:radiation), model.parameterizations)
    dry_model = PrimitiveDryModel(spectral_grid)
    @test findfirst(==(:boundary_layer), dry_model.parameterizations) < findfirst(==(:radiation), dry_model.parameterizations)
end
