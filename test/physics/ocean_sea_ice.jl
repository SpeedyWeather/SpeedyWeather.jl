@testset "Ocean and sea ice models" begin
    
    spectral_grid = SpectralGrid(trunc=31, nlayers=5)

    # just test that these parameters can be set
    SlabOcean(spectral_grid)
    ThermodynamicSeaIce(spectral_grid, temp_freeze=-1.8)

    @testset for OceanModel in ( SeasonalOceanClimatology, 
                            ConstantOceanClimatology,
                            AquaPlanet,
                            SlabOcean)

        @testset for SeaIceModel in (ThermodynamicSeaIce,
                                        Nothing)

            ocean = OceanModel(spectral_grid)
            sea_ice = SeaIceModel(spectral_grid)
            albedo = OceanLandAlbedo(ocean=OceanSeaIceAlbedo(spectral_grid), land=AlbedoClimatology(spectral_grid))

            model = PrimitiveDryModel(spectral_grid; ocean, sea_ice, albedo)
            model.feedback.verbose = false
            simulation = initialize!(model, time=DateTime(2000, 5, 1))
            run!(simulation, period=Day(3))

            @test simulation.model.feedback.nans_detected == false


            @test haskey(simulation.prognostic_variables.ocean, :sea_surface_temperature)

            # Some SSTs may contain NaNs
            # @test all(0 .<= simulation.prognostic_variables.ocean.sea_surface_temperature .<= 330)
            sst = simulation.prognostic_variables.ocean.sea_surface_temperature
            for ij in eachindex(sst)
                if !isnan(sst[ij])
                    @test 0 <= sst[ij] <= 330
                end
            end

            if sea_ice isa Nothing
                @test !haskey(simulation.prognostic_variables.ocean, :sea_ice_concentration)
            else
                @test all(0 .<= simulation.prognostic_variables.ocean.sea_ice_concentration .<= 1)
            end 
        end
    end
end