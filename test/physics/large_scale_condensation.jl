@testset "Large-scale condensation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for r in (0.9, 1.0)
        @testset for snow in (true, false)
            @testset for reevaporation in (0, 10, 30)
                large_scale_condensation = ImplicitCondensation(spectral_grid;
                    relative_humidity_threshold=r, snow, reevaporation)
                model = PrimitiveWetModel(spectral_grid; large_scale_condensation)
                simulation = initialize!(model)
                run!(simulation, period=Day(5))
                
                precip = simulation.diagnostic_variables.physics.precip_large_scale
                snow_fall = simulation.diagnostic_variables.physics.snow_large_scale

                for ij in eachindex(precip)
                    @test precip[ij] >= 0   # precipitation should always be non-negative
                end
            
                if snow
                    for ij in eachindex(snow_fall)
                        @test snow_fall[ij] >= 0     # snow should always be non-negative
                    end

                    @test any(snow_fall .> 0)        # should have some snow
                    @test any(snow_fall .== 0)       # should have some areas without snow

                else
                    @test all(snow_fall .== 0)       # no snow should be produced
                end
            end
        end
    end
end