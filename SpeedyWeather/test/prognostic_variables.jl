@testset "Variables initialize, zero, fill!, one" begin

    NF = Float32
    nlayers = 2
    spectral_grid = SpectralGrid(; NF, nlayers, Grid = FullGaussianGrid)
    model = PrimitiveWetModel(spectral_grid)
    add!(model, Tracer(:abc))
    simulation = initialize!(model)

    # evolve a bit to have nonzero elements
    run!(simulation, period = Day(1), output = false)

    # zero test
    progn = simulation.variables.prognostic

    progn_new = zero(progn)

    for i in eachindex(progn_new.vorticity)
        @test all(progn_new.vorticity[i] .== zero(NF))
        @test all(progn_new.divergence[i] .== zero(NF))
        @test all(progn_new.temperature[i] .== zero(NF))
        @test all(progn_new.humidity[i] .== zero(NF))
        @test all(progn_new.pressure[i] .== zero(NF))
    end

    @test keys(progn_new.tracers) == keys(progn.tracers)

    # copy test
    copy!(progn_new, progn)

    # NaN != NaN, and there's NaNs in the ocean and land, that's why we can't directly test for equality
    # we don't test really all fields, it would just repeat the code that's already there in the main part
    for i in eachindex(progn_new.vorticity)
        @test all(progn_new.vorticity[i] .== progn.vorticity[i])
        @test all(progn_new.divergence[i] .== progn.divergence[i])
        @test all(progn_new.temperature[i] .== progn.temperature[i])
        @test all(progn_new.humidity[i] .== progn.humidity[i])
        @test all(progn_new.pressure[i] .== progn.pressure[i])
    end

    # fill! test
    fill!(progn_new, 1)

    for i in eachindex(progn_new.vorticity)
        @test all(progn_new.vorticity[i] .== one(NF))
        @test all(progn_new.divergence[i] .== one(NF))
        @test all(progn_new.temperature[i] .== one(NF))
        @test all(progn_new.humidity[i] .== one(NF))
        @test all(progn_new.pressure[i] .== one(NF))
    end

    for (key, value) in progn_new.tracers
        for value_i in value
            @test all(value_i .== one(NF))
        end
    end
end
