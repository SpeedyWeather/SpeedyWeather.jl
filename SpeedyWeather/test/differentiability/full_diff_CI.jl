# Enzyme and Julia 1.11 still has some problems, and the test below is broken
# in Julia 1.11
using Enzyme, FiniteDifferences

if VERSION <= v"1.11.0"
    @testset "Complete Differentiability" begin
        # We do extensive correctness checks and tests on the differentiability
        # in a seperate test set. But we do want to ensure in the regular CI that
        # we don't commit some kind of problem for the Enzyme differentiability
        # so, we test here if we get a non-zero gradient from the timestepping.
        spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)          # define resolution
        model = PrimitiveWetModel(; spectral_grid)   # construct model
        simulation = initialize!(model)
        initialize!(simulation)
        run!(simulation, period = Hour(6))

        (; variables, model) = simulation
        (; Δt, Δt_millisec) = model.time_stepping
        dt = 2Δt

        vars = variables
        vars_copy = deepcopy(vars)

        d_vars = make_zero(vars)
        d_model = make_zero(model)

        vars_new = make_zero(vars)
        dvars_new = make_zero(vars)

        # seed dvars_new with ones (output seed)
        for k in keys(dvars_new.prognostic)
            field = getfield(dvars_new.prognostic, k)
            if field isa AbstractArray
                field .= one(eltype(field))
            end
        end

        autodiff(Reverse, timestep_oop!, Const, Duplicated(vars_new, dvars_new), Duplicated(vars, d_vars), Const(dt), Duplicated(model, d_model))
        @test sum(to_vec(d_vars)[1]) != 0

    end
else
    @testset "Complete Differentiability" begin
        @test_broken false
    end
end
