@testset "GPU tracers: fused (default) as one by one" begin
    arch = SpeedyWeather.GPU()

    # 3 tracers fit into one transform, 4 and 11 not (into the largest transform batch of 10 and 19 layers)
    for (Model, nlayers, ntracers) in (
            (ShallowWaterModel, 1, 3), (ShallowWaterModel, 1, 4), (PrimitiveWetModel, 8, 3), (PrimitiveDryModel, 2, 11),
        )
        spectral_grid = SpectralGrid(truncation = 32; nlayers, architecture = arch)
        names = ntuple(i -> Symbol(:tracer, i), ntracers)

        simulations = map((nothing, false)) do fuse_tracers     # default (= fused on GPU) and one by one
            tracers = NamedTuple{names}(ntuple(_ -> Tracer(), ntracers))
            model = isnothing(fuse_tracers) ? Model(spectral_grid; tracers) : Model(spectral_grid; tracers, fuse_tracers)
            model.feedback.verbose = false
            return initialize!(model)
        end
        fused, per_tracer = simulations
        @test SpeedyWeather.fused_tracers(fused.variables, fused.model)
        @test !SpeedyWeather.fused_tracers(per_tracer.variables, per_tracer.model)

        set!(
            fused; namespace = :tracers,
            (name => ((λ, φ, σ) -> exp(-(λ - 36i)^2 / 30^2) * (1 + i * σ) + cosd(φ)) for (i, name) in enumerate(names))...
        )
        copy!(per_tracer.variables, fused.variables)
        tracer0 = Array(fused.variables.prognostic.tracers.tracer1.data)

        run!(fused, steps = 10)
        run!(per_tracer, steps = 10)
        @test fused.model.feedback.nans_detected == false

        # the same up to rounding errors (transforms of other batch sizes)
        for name in names, group in (:prognostic, :grid)
            tracer_fused = Array(getfield(getfield(fused.variables, group).tracers, name).data)
            tracer_per_tracer = Array(getfield(getfield(per_tracer.variables, group).tracers, name).data)
            @test tracer_fused ≈ tracer_per_tracer rtol = 1.0e-4
        end
        @test Array(fused.variables.prognostic.tracers.tracer1.data) != tracer0      # advected
    end
end
