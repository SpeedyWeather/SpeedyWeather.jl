# The Radiation bundle must reproduce the two separate shortwave/longwave components bitwise.

# recursive == that also looks inside Ref (ScalarDim variables) and NamedTuple namespaces
same(a::Base.RefValue, b::Base.RefValue) = a[] == b[]
same(a::NamedTuple, b::NamedTuple) = keys(a) == keys(b) && all(same(a[k], b[k]) for k in keys(a))
same(a, b) = a == b

function init_radiation_state!(vars, model)
    vars.grid.temperature .= 280
    haskey(vars.grid, :humidity) && (vars.grid.humidity .= 1.0e-3)
    vars.grid.pressure .= 1.0e5
    vars.parameterizations.surface_pressure .= 1.0e5
    haskey(vars.parameterizations, :cloud_top) && (vars.parameterizations.cloud_top .= model.spectral_grid.nlayers + 1)
    haskey(vars.parameterizations, :rain_rate) && (vars.parameterizations.rain_rate .= 0)
    vars.parameterizations.ocean.albedo .= 0.5
    vars.parameterizations.land.albedo .= 0.3
    vars.prognostic.ocean.sea_surface_temperature .= 290
    vars.prognostic.land.soil_temperature .= 285
    SpeedyWeather.parameterization!(vars, model.solar_zenith, model)
    return nothing
end

@testset "Radiation bundle" begin
    shortwaves = (TransparentShortwave, OneBandShortwave, OneBandGreyShortwave, Nothing)
    longwaves = (UniformCooling, JeevanjeeRadiation, OneBandLongwave, OneBandGreyLongwave, Nothing)

    @testset "Bit-identical to separate schemes" begin
        # full matrix in Float32 (the default), spot check of the default pair in Float64
        cases = [(Float32, SW, LW) for SW in shortwaves for LW in longwaves]
        push!(cases, (Float64, OneBandShortwave, OneBandLongwave))

        @testset for (NF, SW, LW) in cases
            spectral_grid = SpectralGrid(; NF, truncation = 32, nlayers = 8)
            radiation = Radiation(spectral_grid; shortwave = SW(spectral_grid), longwave = LW(spectral_grid))
            model = PrimitiveWetModel(spectral_grid; radiation)
            initialize!(model.radiation, model)

            vars_separate = Variables(model)
            init_radiation_state!(vars_separate, model)
            vars_bundled = deepcopy(vars_separate)

            for ij in 1:model.spectral_grid.npoints
                SpeedyWeather.parameterization!(ij, vars_separate, radiation.shortwave, model)
                SpeedyWeather.parameterization!(ij, vars_separate, radiation.longwave, model)
            end

            for ij in 1:model.spectral_grid.npoints
                SpeedyWeather.parameterization!(ij, vars_bundled, radiation, model)
            end

            @test vars_separate.tendencies.grid.temperature == vars_bundled.tendencies.grid.temperature
            @test same(vars_separate.parameterizations, vars_bundled.parameterizations)

            # and the tendency is non-trivial unless both streams are off
            if !(SW === Nothing && LW === Nothing)
                @test any(!=(0), vars_bundled.tendencies.grid.temperature)
            end
        end
    end

    @testset "Variables are the union of both streams" begin
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
        model = PrimitiveWetModel(spectral_grid)
        ids(vars) = Set(SpeedyWeather.identifier(v) for v in vars)

        separate = ids(
            (
                SpeedyWeather.variables(model.radiation.shortwave, model)...,
                SpeedyWeather.variables(model.radiation.longwave, model)...,
            )
        )
        @test ids(SpeedyWeather.variables(model.radiation, model)) == separate

        # nothing for one stream still allocates the diagnostics of the other
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; shortwave = nothing))
        vars = Variables(model)
        @test haskey(vars.parameterizations, :outgoing_longwave)
        @test !haskey(vars.parameterizations, :outgoing_shortwave)
    end

    @testset "Model construction and deprecations" begin
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

        # defaults
        model = PrimitiveWetModel(spectral_grid)
        @test model.radiation isa Radiation{<:OneBandShortwave, <:OneBandLongwave}
        @test :radiation in model.parameterizations
        model = PrimitiveDryModel(spectral_grid)
        @test model.radiation.shortwave isa OneBandShortwave     # the grey variant is a OneBandShortwave
        @test model.radiation.longwave isa OneBandLongwave

        # deprecated keywords are wrapped into Radiation
        shortwave = TransparentShortwave(spectral_grid)
        longwave = UniformCooling(spectral_grid)
        model = @test_deprecated PrimitiveWetModel(spectral_grid; shortwave_radiation = shortwave)
        @test model.radiation.shortwave === shortwave
        @test model.radiation.longwave isa OneBandLongwave
        model = @test_deprecated PrimitiveDryModel(spectral_grid; longwave_radiation = longwave)
        @test model.radiation.longwave === longwave
        @test model.radiation.shortwave isa OneBandShortwave

        # deprecated symbols in the parameterizations tuple become a single :radiation
        model = @test_deprecated PrimitiveWetModel(
            spectral_grid;
            parameterizations = (:convection, :shortwave_radiation, :longwave_radiation, :boundary_layer)
        )
        @test model.parameterizations == (:convection, :radiation, :boundary_layer)

        # both old and new keywords together is an error
        @test_throws ArgumentError PrimitiveWetModel(
            spectral_grid;
            radiation = Radiation(spectral_grid), longwave_radiation = longwave
        )

        # initialize! and a few time steps run through
        simulation = initialize!(PrimitiveWetModel(spectral_grid))
        run!(simulation, steps = 2)
        @test all(isfinite, simulation.variables.parameterizations.outgoing_longwave)
    end
end
