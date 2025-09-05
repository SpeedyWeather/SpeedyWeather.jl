@testset "PrognosticVariables initialize, zero, fill!, one" begin
    
    NF = Float32
    nlayers = 2 
    spectral_grid = SpectralGrid(; NF, nlayers, Grid=FullGaussianGrid)
    model = PrimitiveWetModel(spectral_grid)
    add!(model, Tracer(:abc))
    simulation = initialize!(model)
    
    # evolve a bit to have nonzero elements 
    run!(simulation, period=Day(1), output=false)

    # zero test 
    progn = simulation.prognostic_variables

    progn_new = zero(progn)

    for i in eachindex(progn_new.vor)
        @test all(progn_new.vor[i] .== zero(NF))
        @test all(progn_new.div[i] .== zero(NF))
        @test all(progn_new.temp[i] .== zero(NF))
        @test all(progn_new.humid[i] .== zero(NF))
        @test all(progn_new.pres[i] .== zero(NF))
    end

    @test keys(progn_new.tracers) == keys(progn.tracers)

    # copy test 
    copy!(progn_new, progn)

    # NaN != NaN, and there's NaNs in the ocean and land, that's why we can't directly test for equality 
    # we don't test really all fields, it would just repeat the code that's already there in the main part
    for i in eachindex(progn_new.vor)
        @test all(progn_new.vor[i] .== progn.vor[i])
        @test all(progn_new.div[i] .== progn.div[i])
        @test all(progn_new.temp[i] .== progn.temp[i])
        @test all(progn_new.humid[i] .== progn.humid[i])
        @test all(progn_new.pres[i] .== progn.pres[i])
    end

    # fill! test 
    fill!(progn_new, 1)

    for i in eachindex(progn_new.vor)
        @test all(progn_new.vor[i] .== one(NF))
        @test all(progn_new.div[i] .== one(NF))
        @test all(progn_new.temp[i] .== one(NF))
        @test all(progn_new.humid[i] .== one(NF))
        @test all(progn_new.pres[i] .== one(NF))
    end
    
    for (key, value) in progn_new.tracers
        for value_i in value 
            @test all(value_i .== one(NF))
        end 
    end 
end