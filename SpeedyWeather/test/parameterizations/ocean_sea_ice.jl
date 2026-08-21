@testset "Ocean and sea ice models" begin

    spectral_grid = SpectralGrid(truncation = 32, nlayers = 5)

    # just test that these parameters can be set
    SlabOcean(spectral_grid)
    ThermodynamicSeaIce(spectral_grid, freezing_temperature = 273.15)

    @testset for OceanModel in (
            SeasonalOceanClimatology,
            ConstantOceanClimatology,
            AquaPlanet,
            SlabOcean,
            PrescribedOcean,
        )

        @testset for SeaIceModel in (
                ThermodynamicSeaIce,
                PrescribedSeaIce,
                Nothing,
            )

            ocean = OceanModel(spectral_grid)
            sea_ice = SeaIceModel(spectral_grid)
            albedo = OceanLandAlbedo(ocean = OceanSeaIceAlbedo(spectral_grid), land = AlbedoClimatology(spectral_grid))

            model = PrimitiveDryModel(spectral_grid; ocean, sea_ice, albedo)
            model.feedback.verbose = false
            simulation = initialize!(model, time = DateTime(2000, 5, 1))

            # PrescribedOcean does not initialize SST set to something reasonable if it is all zeros
            if all(simulation.variables.prognostic.ocean.sea_surface_temperature .== 0)
                simulation.variables.prognostic.ocean.sea_surface_temperature .= 285
            end

            run!(simulation, period = Day(1))

            @test simulation.model.feedback.nans_detected == false
            @test haskey(simulation.variables.prognostic.ocean, :sea_surface_temperature)

            @test ArrayDimensions.hastime(simulation.variables.prognostic.ocean.sea_surface_temperature) == (ocean isa SpeedyWeather.AbstractDynamicOcean)
            sst = get_step(simulation.variables.prognostic.ocean.sea_surface_temperature, 1)
            @test ndims(sst) == 1
            @test all(0 .<= sst .<= 330)

            if sea_ice isa Nothing
                @test !haskey(simulation.variables.prognostic.ocean, :sea_ice_concentration)
            else
                @test ArrayDimensions.hastime(simulation.variables.prognostic.ocean.sea_ice_concentration) == (sea_ice isa SpeedyWeather.AbstractDynamicSeaIce)
                sic = get_step(simulation.variables.prognostic.ocean.sea_ice_concentration, 1)
                @test ndims(sic) == 1
                @test all(0 .<= sic .<= 1)
            end
        end
    end
end
