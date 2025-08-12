@testset "Ocean and sea ice models" begin
    
    spectral_grid = SpectralGrid(trunc=31, nlayers=5)

    for ocean in (  SeasonalOceanClimatology(spectral_grid), 
                    ConstantOceanClimatology(spectral_grid),
                    AquaPlanet(spectral_grid),
                    SlabOcean(spectral_grid, sea_ice_insulation = (x) -> x))

        for sea_ice in (ThermodynamicSeaIce(spectral_grid, temp_freeze=0),
                        NoSeaIce(spectral_grid))

            albedo = Albedo(ocean=OceanSeaIceAlbedo(spectral_grid, land=AlbedoClimatology(spectral_grid)))

            model = PrimitiveDryModel(spectral_grid; ocean, sea_ice, albedo)
            model.feedback.verbose = false
            simulation = initialize!(model, time=DateTime(2000,5,1))
            run!(simulation, period=Day(3))

            @test simulation.model.feedback.nans_detected == false

            @test all(0 .<= simulation.prognostic_variables.ocean.sea_ice_concentration .<= 1)
            @test all(0 .<= simulation.diagnostic_variables.physics.ocean.albedo .<= 1)

            if sea_ice isa NoSeaIce
                @test all(simulation.prognostic_variables.ocean.sea_ice_concentration .== 0)
            end
        end
    end
end