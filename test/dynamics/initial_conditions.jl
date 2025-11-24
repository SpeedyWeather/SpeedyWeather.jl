@testset "initial_conditions" begin

    @testset "ZeroInitially" begin
        spectral_grid = SpectralGrid(nlayers = 1)
        ic = ZeroInitially()
        model = BarotropicModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)

        # Check that vorticity is initialized to zero
        vor = prognostic_variables.vor
        @test all(vor .== 0)
        @test !any(isnan.(vor))
    end

    @testset "RandomVorticity" begin
        spectral_grid = SpectralGrid(nlayers = 1)

        # Test with default parameters
        ic = RandomVorticity(spectral_grid)
        model = BarotropicModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        vor = prognostic_variables.vor

        # Check vorticity is non-zero and has no NaNs
        @test !all(vor .== 0)
        @test !any(isnan.(vor))

        # Test different parameters produce different results
        ic_custom = RandomVorticity(spectral_grid; power = -2, amplitude = 1.0e-5, max_wavenumber = 15)
        model_custom = BarotropicModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        vor_custom = prognostic_variables_custom.vor

        @test !any(isnan.(vor_custom))
        @test vor_custom != vor  # Different parameters should give different results
    end

    @testset "RandomVelocity" begin
        spectral_grid = SpectralGrid(nlayers = 1)

        # Test with default parameters
        ic = RandomVelocity(spectral_grid)
        model = BarotropicModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        vor = prognostic_variables.vor

        # Check vorticity is non-zero and has no NaNs
        @test !all(vor .== 0)
        @test !any(isnan.(vor))

        # Test different parameters produce different results
        ic_custom = RandomVelocity(spectral_grid; max_speed = 40, truncation = 10)
        model_custom = BarotropicModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        vor_custom = prognostic_variables_custom.vor

        @test !any(isnan.(vor_custom))
        @test vor_custom != vor  # Different parameters should give different results
    end

    @testset "ZonalJet" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with default parameters (Galewsky jet)
        ic = ZonalJet(spectral_grid)
        model = ShallowWaterModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        vor = prognostic_variables.vor
        pres = prognostic_variables.pres

        # Check vorticity and pressure are initialized
        @test !all(vor .== 0)
        @test !any(isnan.(vor))
        @test !all(pres .== 0)
        @test !any(isnan.(pres))

        # Test different parameters produce different results
        ic_custom = ZonalJet(spectral_grid; latitude = 30, width = 15, umax = 60)
        model_custom = ShallowWaterModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        vor_custom = prognostic_variables_custom.vor

        @test !any(isnan.(vor_custom))
        @test vor_custom != vor  # Different parameters should give different results
    end

    @testset "ZonalWind" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with default parameters (Jablonowski)
        ic = ZonalWind(spectral_grid)
        model = PrimitiveDryModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        vor = prognostic_variables.vor
        div = prognostic_variables.div

        # Check vorticity and divergence are initialized
        @test !all(vor .== 0)
        @test !any(isnan.(vor))
        @test !all(div .== 0)
        @test !any(isnan.(div))

        # Test different parameters produce different results
        ic_custom = ZonalWind(spectral_grid; u₀ = 40, perturb_uₚ = 2)
        model_custom = PrimitiveDryModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        vor_custom = prognostic_variables_custom.vor

        @test !any(isnan.(vor_custom))
        @test vor_custom != vor  # Different parameters should give different results
    end

    @testset "RossbyHaurwitzWave" begin
        spectral_grid = SpectralGrid(nlayers = 1)

        # Test with ShallowWater model
        ic = RossbyHaurwitzWave(spectral_grid)
        model = ShallowWaterModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        vor = prognostic_variables.vor
        pres = prognostic_variables.pres

        # Check vorticity and pressure are initialized
        @test !all(vor .== 0)
        @test !any(isnan.(vor))
        @test !all(pres .== 0)
        @test !any(isnan.(pres))

        # Test different parameters produce different results
        ic_custom = RossbyHaurwitzWave(spectral_grid; m = 6)
        model_custom = ShallowWaterModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        vor_custom = prognostic_variables_custom.vor

        @test !any(isnan.(vor_custom))
        @test vor_custom != vor  # Different wavenumber should give different results
    end

    @testset "JablonowskiTemperature" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with default parameters
        ic = JablonowskiTemperature(spectral_grid)
        model = PrimitiveDryModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        temp = prognostic_variables.temp

        # Check temperature is initialized
        @test !all(temp .== 0)
        @test !any(isnan.(temp))

        # Test different parameters produce different results
        ic_custom = JablonowskiTemperature(spectral_grid; ΔT = 4.8e5)
        model_custom = PrimitiveDryModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        temp_custom = prognostic_variables_custom.temp

        @test !any(isnan.(temp_custom))
        @test temp_custom != temp  # Different parameters should give different results
    end

    @testset "PressureOnOrography" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with PrimitiveDry model
        ic = PressureOnOrography(spectral_grid)
        model = PrimitiveDryModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        pres = prognostic_variables.pres

        # Check pressure is initialized
        @test !all(pres .== 0)
        @test !any(isnan.(pres))

        # Test with PrimitiveWet model
        model_wet = PrimitiveWetModel(spectral_grid)
        prognostic_variables_wet = PrognosticVariables(model_wet)
        initialize!(prognostic_variables_wet, ic, model_wet)
        pres_wet = prognostic_variables_wet.pres

        @test !all(pres_wet .== 0)
        @test !any(isnan.(pres_wet))

        @test pres_wet != pres  # Different models should give different results
    end

    @testset "ConstantPressure" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with PrimitiveDry model
        ic = ConstantPressure(spectral_grid)
        model = PrimitiveDryModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        pres = prognostic_variables.pres

        # Check pressure is initialized to constant value
        @test !any(isnan.(pres))

        # Test with ShallowWater model (should do nothing)
        spectral_grid_sw = SpectralGrid(nlayers = 1)
        ic_sw = ConstantPressure(spectral_grid_sw)
        model_sw = ShallowWaterModel(spectral_grid_sw)
        prognostic_variables_sw = PrognosticVariables(model_sw)
        initialize!(prognostic_variables_sw, ic_sw, model_sw)
        pres_sw = prognostic_variables_sw.pres

        @test all(pres_sw .≈ 0)
        @test !any(isnan.(pres_sw))
    end

    @testset "ConstantRelativeHumidity" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with default relative humidity
        ic = ConstantRelativeHumidity(spectral_grid)
        model = PrimitiveWetModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        humid = prognostic_variables.humid

        # Check humidity is initialized
        @test !all(humid .== 0)
        @test !any(isnan.(humid))

        # Test different parameters produce different results
        ic_custom = ConstantRelativeHumidity(spectral_grid; relhumid_ref = 0.5)
        model_custom = PrimitiveWetModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        humid_custom = prognostic_variables_custom.humid

        @test !any(isnan.(humid_custom))
        @test humid_custom != humid  # Different relative humidity should give different results
    end

    @testset "RandomWaves" begin
        spectral_grid = SpectralGrid(nlayers = 1)

        # Test with default parameters
        ic = RandomWaves(spectral_grid)
        model = ShallowWaterModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        pres = prognostic_variables.pres

        # Check pressure is initialized
        @test !all(pres .== 0)
        @test !any(isnan.(pres))

        # Test different parameters produce different results
        ic_custom = RandomWaves(spectral_grid; A = 1000, lmin = 5, lmax = 20)
        model_custom = ShallowWaterModel(spectral_grid)
        prognostic_variables_custom = PrognosticVariables(model_custom)
        initialize!(prognostic_variables_custom, ic_custom, model_custom)
        pres_custom = prognostic_variables_custom.pres

        @test !any(isnan.(pres_custom))
        @test pres_custom != pres  # Different parameters should give different results
    end

    @testset "StartFromRest" begin
        spectral_grid = SpectralGrid(nlayers = 8)

        # Test with PrimitiveDry model
        ic = StartFromRest(spectral_grid)
        model = PrimitiveDryModel(spectral_grid)
        prognostic_variables = PrognosticVariables(model)
        initialize!(prognostic_variables, ic, model)
        pres = prognostic_variables.pres
        temp = prognostic_variables.temp

        # Check pressure and temperature are initialized
        @test !any(isnan.(pres))
        @test !any(isnan.(temp))

        # Test with PrimitiveWet model
        model_wet = PrimitiveWetModel(spectral_grid)
        prognostic_variables_wet = PrognosticVariables(model_wet)
        initialize!(prognostic_variables_wet, ic, model_wet)
        humid = prognostic_variables_wet.humid

        @test !any(isnan.(humid))
    end
end
