# Tests for the CUDA-Graphs accelerated batched Fourier transform (SpeedyTransformsCUDAExt).
# Runs a GPU model and checks the graph cache stays bounded: the time loop reuses the same
# variable buffers every step, so each (size, field) pair is captured once and replayed —
# the per-direction exec dicts must stay well under MAX_GRAPHS and not grow with more steps.
# Only executes when the CUDA extension is actually loaded.

@testset "CUDA Graphs: bounded graph cache over a GPU model run" begin
    ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsCUDAExt)
    if ext !== nothing
        # CUDA-Graphs path is on by default (SpectralTransform's cuda_graphs = true)
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

# Focused regression test for the graph-cache keying. The time stepping fetches each transform's
# grid field via a per-step view (`get_prognostic_step`/`get_tendency_step` → `get_step` →
# `field_view`). On GPU every such view is a FRESH `CuArray` wrapper aliasing the same device
# memory, so the cache must key on the stable device pointer — not the churning wrapper identity —
# or it captures a new graph every step and grows without bound (the failure caught by the
# full-model test above, isolated here without a model run).
@testset "CUDA Graphs: per-step views of one buffer reuse a single graph" begin
    ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsCUDAExt)
    if ext !== nothing
        spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = GPU())
        S = SpectralTransform(spectral_grid)
        nlayers = spectral_grid.nlayers

        # grid field WITH a leapfrog-style step dimension (npoints, nlayers, nsteps), like the new
        # time stepping's prognostic grid variables; each get_step is a fresh view of step 2
        gridded = rand(Float32, spectral_grid.grid, nlayers, 2)
        specs = rand(ComplexF32, spectral_grid.spectrum, nlayers)

        # the hazard: the per-step view wrapper identity is NOT stable across calls
        @test get_step(gridded, 2).data !== get_step(gridded, 2).data

        # inverse (spectral→grid) into a fresh step-2 view each call → only ONE graph captured
        ext.clear_fourier_graph_cache!()
        for _ in 1:4
            transform!(get_step(gridded, 2), specs, S)
        end
        @test sum(length(c.inverse_execs) for c in values(ext.GRAPH_CACHES); init = 0) == 1

        # forward (grid→spectral) from a fresh step-2 view each call → only ONE graph captured
        ext.clear_fourier_graph_cache!()
        for _ in 1:4
            transform!(specs, get_step(gridded, 2), S)
        end
        @test sum(length(c.forward_execs) for c in values(ext.GRAPH_CACHES); init = 0) == 1

        ext.clear_fourier_graph_cache!()
    end
end
