@testset "copy!(::Variables, ::Variables)" begin

    @testset "Barotropic" begin
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
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
        @test vars_copy.prognostic.clock.step_counter == vars.prognostic.clock.step_counter
        @test vars_copy.prognostic.clock.time_step_counter == vars.prognostic.clock.time_step_counter
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
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
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
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 5)
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
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
        model = BarotropicModel(spectral_grid)
        model.feedback.verbose = false
        simulation = initialize!(model)

        vars = simulation.variables
        vars_copy = Variables(model)
        result = copy!(vars_copy, vars)
        @test result === vars_copy
    end
end

# The groups `copy!` copies (not the scratch) on the CPU, to compare with `==`. Without the fused parents
# (compared through their views) and the clock and scale, which `==` compares by identity (mutable).
# Scalar variables (e.g. greenhouse gases) are `Ref`s, also compared by identity, so unwrap them to their value
unwrap_refs(x) = x
unwrap_refs(x::Base.RefValue) = x[]
unwrap_refs(x::NamedTuple) = map(unwrap_refs, x)
variables_on_cpu(vars) = map((:prognostic, :grid, :tendencies, :dynamics, :parameterizations, :particles)) do group
    unwrap_refs(on_architecture(SpeedyWeather.CPU(), Base.structdiff(getfield(vars, group), NamedTuple{(:clock, :scale)})))
end

# On JLArrays (as on CUDA) a contiguous view of a fused parent is not a `SubArray` but another JLArray sharing
# the parent's memory: `copy!` has to skip it (copying its parent instead), `copyto!` between it and a `SubArray`
# on the CPU would fall back to scalar indexing
using JLArrays

@testset "copy! between CPU and JLArray variables" begin
    spectral_grid = SpectralGrid(truncation = 11, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps = 3)                     # non-trivial data in all groups
    vars = simulation.variables
    reference = variables_on_cpu(vars)

    architecture = SpeedyWeather.GPU(JLArrays.JLBackend())
    jlarray_model = PrimitiveWetModel(SpectralGrid(truncation = 11, nlayers = 4; architecture))
    vars_jlarray = copy!(Variables(jlarray_model), vars)        # CPU → JLArray
    @test variables_on_cpu(vars_jlarray) == reference
    vars_jlarray2 = copy!(Variables(jlarray_model), vars_jlarray)   # JLArray → JLArray
    @test variables_on_cpu(vars_jlarray2) == reference
    vars_back = copy!(Variables(model), vars_jlarray2)          # JLArray → CPU
    @test variables_on_cpu(vars_back) == reference
end
