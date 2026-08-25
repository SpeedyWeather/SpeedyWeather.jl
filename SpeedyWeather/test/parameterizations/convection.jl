@testset "Convection" begin

    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    @testset for Convection in (
            Nothing,
            BettsMillerConvection,
            BettsMillerDryConvection,
            ConvectiveHeating,
        )

        @testset for Model in (
                PrimitiveDryModel,
                PrimitiveWetModel,
            )

            # that combination is not defined
            if ~(Convection == BettsMillerConvection && Model == PrimitiveDryModel)
                convection = Convection(spectral_grid)
                model = Model(spectral_grid; convection)
                model.feedback.verbose = false
                simulation = initialize!(model)

                run!(simulation, steps = 36)
            end
        end
    end

    @testset "Entrainment functors" begin
        for NF in (Float32, Float64)
            ne = NoEntrainment()
            for σ in NF.((0, 0.3, 0.7, 1))
                @test ne(σ) == 0
                @test (@inferred ne(σ)) == 0
            end

            LinearEntrainment(NF)                   # constructor
            le = LinearEntrainment{NF}(; σ_entrainment = 0.5, surface_entrainment = 0.8)
            @test le(NF(1)) ≈ NF(0.8)               # surface: full surface_entrainment
            @test le(NF(0.5)) == 0                  # at σ_entrainment: zero
            @test le(NF(0.2)) == 0                  # above σ_entrainment: zero
            @test le(NF(0.75)) < le(NF(1))          # monotonically non-decreasing in σ
            @test le(NF(0.5)) <= le(NF(0.75)) <= le(NF(1))
            @test (@inferred le(NF(0.8))) isa NF

            ConstantEntrainment(NF)                 # constructor only
            ce = ConstantEntrainment{NF}(; entrainment_rate = 0.3)
            for σ in NF.((0, 0.3, 0.7, 1))
                @test ce(σ) ≈ NF(0.3)
            end
            @test (@inferred ce(NF(0.5))) isa NF
        end
    end

    @testset "Entrainment in a model" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        for entrainment in (NoEntrainment(), LinearEntrainment(spectral_grid), ConstantEntrainment(spectral_grid))
            wet_model = PrimitiveWetModel(spectral_grid; convection = BettsMillerConvection(spectral_grid; entrainment))
            wet_model.feedback.verbose = false
            wet_simulation = initialize!(wet_model)
            run!(wet_simulation, steps = 12)
            @test !any(isnan, wet_simulation.variables.prognostic.temperature)

            dry_model = PrimitiveDryModel(spectral_grid; convection = BettsMillerDryConvection(spectral_grid; entrainment))
            dry_model.feedback.verbose = false
            dry_simulation = initialize!(dry_model)
            run!(dry_simulation, steps = 12)
            @test !any(isnan, dry_simulation.variables.prognostic.temperature)
        end
    end

    @testset "Entrainment has an effect" begin
        # direct, deterministic comparison on the *same* atmospheric snapshot (rather than
        # comparing two independent model runs, which diverge chaotically regardless of
        # entrainment and so cannot isolate its effect)
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        model = PrimitiveWetModel(spectral_grid)
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 5)

        vars = simulation.variables
        nlayers = spectral_grid.nlayers
        temp = SpeedyWeather.get_prognostic_step(vars.grid.temperature, model.time_stepping, model.convection)
        humid = SpeedyWeather.get_prognostic_step(vars.grid.humidity, model.time_stepping, model.convection)
        geopotential = vars.dynamics.geopotential
        pₛ = vars.parameterizations.surface_pressure
        σ = model.geometry.σ_levels_full
        RH = model.convection.relative_humidity

        temp_ref_none = similar(vars.scratch.grid.a)
        humid_ref_none = similar(vars.scratch.grid.b)
        temp_ref_ent = similar(vars.scratch.grid.a)
        humid_ref_ent = similar(vars.scratch.grid.b)
        entrainment = ConstantEntrainment(spectral_grid; entrainment_rate = 0.3)

        any_different = false
        for ij in 1:spectral_grid.npoints
            lzb_none = SpeedyWeather.pseudo_adiabat!(
                ij, temp_ref_none, humid_ref_none, temp, humid, geopotential,
                pₛ[ij], σ, model.atmosphere, RH, NoEntrainment(),
            )
            lzb_ent = SpeedyWeather.pseudo_adiabat!(
                ij, temp_ref_ent, humid_ref_ent, temp, humid, geopotential,
                pₛ[ij], σ, model.atmosphere, RH, entrainment,
            )
            # entrainment dilutes the rising parcel with (cooler, drier) environmental air,
            # which can only reduce or maintain its buoyancy, never increase the depth of
            # convection: the LZB with entrainment is never above (i.e. never a smaller
            # layer index than) the LZB without entrainment, for the same environment
            @test lzb_ent >= lzb_none
            any_different |= temp_ref_ent[ij, :] != temp_ref_none[ij, :]
        end
        @test any_different   # entrainment must actually change the reference profile somewhere
    end

    @testset "Reference profile invariants" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        model = PrimitiveWetModel(spectral_grid)
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 5)

        vars = simulation.variables
        nlayers = spectral_grid.nlayers
        temp = SpeedyWeather.get_prognostic_step(vars.grid.temperature, model.time_stepping, model.convection)
        humid = SpeedyWeather.get_prognostic_step(vars.grid.humidity, model.time_stepping, model.convection)
        geopotential = vars.dynamics.geopotential
        pₛ = vars.parameterizations.surface_pressure
        σ = model.geometry.σ_levels_full
        temp_ref_profile = vars.scratch.grid.a
        humid_ref_profile = vars.scratch.grid.b

        for ij in (1, spectral_grid.npoints ÷ 2, spectral_grid.npoints)
            lzb = SpeedyWeather.pseudo_adiabat!(
                ij, temp_ref_profile, humid_ref_profile, temp, humid, geopotential,
                pₛ[ij], σ, model.atmosphere, model.convection.relative_humidity, model.convection.entrainment,
            )
            @test 1 <= lzb <= nlayers
            @test !any(isnan, view(temp_ref_profile, ij, :))
            @test !any(isnan, view(humid_ref_profile, ij, :))
            for k in 1:(lzb - 1)
                @test temp_ref_profile[ij, k] == temp[ij, k]
                @test humid_ref_profile[ij, k] == humid[ij, k]
            end
        end

        temp_ref_profile_dry = vars.scratch.grid.a
        for ij in (1, spectral_grid.npoints ÷ 2, spectral_grid.npoints)
            temp_parcel = temp[ij, nlayers]
            lzb = SpeedyWeather.dry_adiabat!(
                ij, temp_ref_profile_dry, temp, temp_parcel, σ, model.atmosphere, NoEntrainment(),
            )
            @test 1 <= lzb <= nlayers
            @test !any(isnan, view(temp_ref_profile_dry, ij, :))
            for k in 1:(lzb - 1)
                @test temp_ref_profile_dry[ij, k] == temp[ij, k]
            end
        end
    end

    @testset "Convective snow" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        model = PrimitiveWetModel(spectral_grid)  # snow = true by default
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 20)   # spin up so some columns actually convect

        vars = simulation.variables
        @test model.convection.snow == true
        @test haskey(vars.parameterizations, :snow_convection)
        @test haskey(vars.parameterizations, :snow_rate_convection)
        @test haskey(vars.parameterizations, :snow_rate)

        # two isolated calls of convection! on the *same* atmospheric state (vars.grid is not
        # mutated by convection!), with snow on vs off, so the comparison isn't confounded by
        # chaotic divergence between separate model runs
        convection_snow = model.convection
        convection_no_snow = BettsMillerConvection(spectral_grid; snow = false)

        vars.parameterizations.rain_convection .= 0
        vars.parameterizations.snow_convection .= 0
        SpeedyWeather._column_parameterizations_cpu!(vars, (convection = convection_snow,), model)
        rain_on = copy(vars.parameterizations.rain_convection)
        snow_on = copy(vars.parameterizations.snow_convection)

        vars.parameterizations.rain_convection .= 0
        vars.parameterizations.snow_convection .= 0
        SpeedyWeather._column_parameterizations_cpu!(vars, (convection = convection_no_snow,), model)
        rain_off = copy(vars.parameterizations.rain_convection)
        snow_off = copy(vars.parameterizations.snow_convection)

        @test all(snow_off .== 0)                    # snow = false: never snows
        @test any(snow_on .> 0)                       # some (polar) column does snow by default
        # snow only swaps rain <-> snow after rain_convection is fully computed, so the total
        # is exactly conserved (not merely approximately) between the two configurations
        @test rain_on .+ snow_on == rain_off
    end
end
