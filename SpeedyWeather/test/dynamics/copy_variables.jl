@testset "copy!(::Variables, ::Variables)" begin

    @testset "Barotropic" begin
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 1)
        model = BarotropicModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Day(1), output = false)

        vars = simulation.variables
        vars_copy = Variables(model)
        copy!(vars_copy, vars)

        # prognostic arrays copied
        @test vars_copy.prognostic.vorticity == vars.prognostic.vorticity
        # clock copied
        @test vars_copy.prognostic.clock.time == vars.prognostic.clock.time
        @test vars_copy.prognostic.clock.timestep_counter == vars.prognostic.clock.timestep_counter
        # scalar Ref copied
        @test vars_copy.prognostic.scale[] == vars.prognostic.scale[]

        # grid variables copied
        @test vars_copy.grid.vorticity == vars.grid.vorticity
        @test vars_copy.grid.u == vars.grid.u
        @test vars_copy.grid.v == vars.grid.v

        # tendencies copied
        @test vars_copy.tendencies.vorticity == vars.tendencies.vorticity

        # verify it's a deep copy, not aliasing
        vars.prognostic.vorticity .= 0
        @test !all(vars_copy.prognostic.vorticity .== 0)
    end

    @testset "ShallowWater" begin
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 1)
        model = ShallowWaterModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Day(1), output = false)

        vars = simulation.variables
        vars_copy = Variables(model)
        copy!(vars_copy, vars)

        @test vars_copy.prognostic.vorticity == vars.prognostic.vorticity
        @test vars_copy.prognostic.divergence == vars.prognostic.divergence
        @test vars_copy.prognostic.η == vars.prognostic.η
        @test vars_copy.grid.vorticity == vars.grid.vorticity
        @test vars_copy.grid.divergence == vars.grid.divergence

        # no aliasing
        vars.prognostic.divergence .= 0
        @test !all(vars_copy.prognostic.divergence .== 0)
    end

    @testset "PrimitiveDry" begin
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 5)
        model = PrimitiveDryModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Day(1), output = false)

        vars = simulation.variables
        vars_copy = Variables(model)
        copy!(vars_copy, vars)

        # prognostic spectral variables
        @test vars_copy.prognostic.vorticity == vars.prognostic.vorticity
        @test vars_copy.prognostic.divergence == vars.prognostic.divergence
        @test vars_copy.prognostic.temperature == vars.prognostic.temperature
        @test vars_copy.prognostic.pressure == vars.prognostic.pressure

        # grid variables
        @test vars_copy.grid.vorticity == vars.grid.vorticity
        @test vars_copy.grid.temperature == vars.grid.temperature

        # namespaced prognostic variables (ocean, land)
        if haskey(vars.prognostic, :ocean)
            for key in keys(vars.prognostic.ocean)
                @test getfield(vars_copy.prognostic.ocean, key) == getfield(vars.prognostic.ocean, key)
            end
        end

        if haskey(vars.prognostic, :land)
            for key in keys(vars.prognostic.land)
                @test getfield(vars_copy.prognostic.land, key) == getfield(vars.prognostic.land, key)
            end
        end

        # dynamics
        @test vars_copy.dynamics.geopotential == vars.dynamics.geopotential
    end

    @testset "copy! returns dest" begin
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 1)
        model = BarotropicModel(spectral_grid)
        simulation = initialize!(model)

        vars = simulation.variables
        vars_copy = Variables(model)
        result = copy!(vars_copy, vars)
        @test result === vars_copy
    end
end
