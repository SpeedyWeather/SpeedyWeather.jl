# Tests for the CUDA-Graphs accelerated batched Fourier transform (SpeedyTransformsCUDAExt).
# Runs a GPU model and checks the graph cache stays bounded: the time loop reuses the same
# variable buffers every step, so each (size, field) pair is captured once and replayed —
# the per-direction exec dicts must stay well under MAX_GRAPHS and not grow with more steps.
# Only executes when the CUDA extension is actually loaded.

@testset "CUDA Graphs: bounded graph cache over a GPU model run" begin
    ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsCUDAExt)
    if ext !== nothing
        ext.FOURIER_GRAPHS_ENABLED[] = true

        spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = GPU())
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)

        # isolate the time loop: drop any graphs captured during initialize!
        ext.clear_fourier_graph_cache!()

        # (max graphs in any one exec dict, total graphs, capture failures) across all caches
        function cache_stats()
            maxlen = total = failed = 0
            for c in values(ext.GRAPH_CACHES), execs in (c.forward_execs, c.inverse_execs)
                maxlen = max(maxlen, length(execs))
                total += length(execs)
                failed += count(e -> e === nothing, values(execs))
            end
            return (; maxlen, total, failed)
        end

        run!(simulation, steps = 5)
        s1 = cache_stats()
        run!(simulation, steps = 5)     # 5 more steps must add no new graphs (buffers reused)
        s2 = cache_stats()

        @test simulation.model.feedback.nans_detected == false
        @test s1.total > 0                  # graphs were actually captured
        @test s2.failed == 0                # nothing fell back to the un-captured direct loop
        @test s1.maxlen < ext.MAX_GRAPHS    # cache stays under the cap
        @test s2.total == s1.total          # extra steps add no graphs → cache is bounded

        ext.clear_fourier_graph_cache!()
    end
end
