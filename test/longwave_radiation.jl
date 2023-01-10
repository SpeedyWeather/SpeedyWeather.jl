@testset "Parametrization: longwave radiation" begin
    @testset "radset!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation)

            # Just check the last band
            @test isapprox(
                model.parameterization_constants.fband[end, :],
                [0.19498351, 0.12541235, 0.33106664, 0.2985375],
                rtol=0.1
            )
        end
    end
    @testset "radlw_down!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev=8)

            nlev = model.parameters.nlev
            nband = model.parameters.nband
            n_stratosphere_levels = model.parameters.n_stratosphere_levels
            column = ColumnVariables{NF}(nlev=nlev, nband=nband, n_stratosphere_levels=n_stratosphere_levels)
            column.temp = fill(300., nlev)
            column.wvi = fill(0.5, nlev, 2)
            column.tau2 = fill(0.5, nlev, 4)
            
            SpeedyWeather.radlw_down!(column, model)

            # Just check fsfcd and dfabs as used by radlw_up! 
            @test column.fsfcd ≈ 447.29436216
            @test column.dfabs ≈ [-71.86727228, -182.21961386, -91.10980693, -45.55490346, -22.77745173, -11.38872587, -5.69436293, -25.35141147]
        end
    end
    @testset "compute_bbe!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev=8)

            nlev = model.parameters.nlev
            nband = model.parameters.nband
            n_stratosphere_levels = model.parameters.n_stratosphere_levels
            column = ColumnVariables{NF}(nlev=nlev, nband=nband, n_stratosphere_levels=n_stratosphere_levels)
            column.ts = 320.
            
            SpeedyWeather.compute_bbe!(column, model)

            @test column.fsfcu ≈ 582.65174016
        end
    end
    @testset "radlw_up!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev=8)

            nlev = model.parameters.nlev
            nband = model.parameters.nband
            n_stratosphere_levels = model.parameters.n_stratosphere_levels
            column = ColumnVariables{NF}(nlev=nlev, nband=nband, n_stratosphere_levels=n_stratosphere_levels)
            column.temp = fill(300., nlev)
            column.wvi = fill(0.5, nlev, 2)
            column.tau2 = fill(0.5, nlev, 4)
            column.ts = 320.
            column.stratc = fill(0.5, 2)
            
            SpeedyWeather.radlw_down!(column, model)
            SpeedyWeather.compute_bbe!(column, model)
            SpeedyWeather.radlw_up!(column, model)

            # Just check what's needed
            @test column.fsfc ≈ 135.357378
            @test column.ftop ≈ 474.76064406
        end
    end
    @testset "longwave_radiation!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev=8)

            nlev = model.parameters.nlev
            nband = model.parameters.nband
            n_stratosphere_levels = model.parameters.n_stratosphere_levels
            column = ColumnVariables{NF}(nlev=nlev, nband=nband, n_stratosphere_levels=n_stratosphere_levels)
            column.temp = fill(300., nlev)
            column.wvi = fill(0.5, nlev, 2)
            column.tau2 = fill(0.5, nlev, 4)
            column.ts = 320.
            column.stratc = fill(0.5, 2)
            
            SpeedyWeather.longwave_radiation!(column, model)

            # Just check what's needed
            @test column.fsfc ≈ 135.357378
            @test column.ftop ≈ 474.76064406
        end
    end
end
