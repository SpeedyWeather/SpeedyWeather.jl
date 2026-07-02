# Shared test logic for CUDA-graphs and HIP-graphs tests.
# Called from cuda_graphs.jl and hip_graphs.jl with the appropriate extension and prefix.
# GRAPH_CACHES, MAX_GRAPHS, and clear_fourier_graph_cache! now live in SpeedyTransforms
# (src/fourier_gpu_graphs.jl); ext is only used to check whether the extension is loaded.
#
# `expect_capture` is false for HIP: AMDGPU's run_graph! never attempts graph capture at all
# (see SpeedyTransformsAMDGPUExt.jl — ROCm's stream-capture validator doesn't reliably reject
# illegal-to-capture operations, so a captured graph can silently replay into invalid memory).
# The allocation-free direct loop still runs on every call, so results must still be correct
# and no graphs should ever appear in the cache.

function test_gpu_graphs(ext, prefix; expect_capture::Bool = true)
    @testset "$prefix Graphs: bounded graph cache over a GPU model run" begin
        if ext !== nothing
            spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = SpeedyWeather.GPU())
            model = PrimitiveWetModel(spectral_grid)
            simulation = initialize!(model)

            SpeedyTransforms.clear_fourier_graph_cache!()

            function cache_stats()
                maxlen = total = failed = 0
                for c in values(SpeedyTransforms.GRAPH_CACHES), execs in (c.forward_execs, c.inverse_execs)
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
            @test s2.failed == 0                # nothing fell back to the un-captured direct loop
            @test s2.total == s1.total          # extra steps add no graphs → cache is bounded
            if expect_capture
                @test s1.total > 0                  # graphs were actually captured
                @test s1.maxlen < SpeedyTransforms.MAX_GRAPHS    # cache stays under the cap
            else
                @test s1.total == 0                 # capture disabled: never captures
            end

            SpeedyTransforms.clear_fourier_graph_cache!()
        end
    end

    # Focused regression test for the graph-cache keying. The time stepping fetches each
    # transform's grid field via a per-step view (`get_step` → `field_view`). On GPU every
    # such view is a FRESH array wrapper aliasing the same device memory, so the cache must
    # key on the stable device pointer — not the churning wrapper identity — or it captures a
    # new graph every step and grows without bound. Not applicable when capture is disabled.
    @testset "$prefix Graphs: per-step views of one buffer reuse a single graph" begin
        if ext !== nothing && expect_capture
            spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = SpeedyWeather.GPU())
            S = SpectralTransform(spectral_grid)
            nlayers = spectral_grid.nlayers

            gridded = rand(Float32, spectral_grid.grid, nlayers, 2)
            specs = rand(ComplexF32, spectral_grid.spectrum, nlayers)

            # the hazard: the per-step view wrapper identity is NOT stable across calls
            @test get_step(gridded, 2).data !== get_step(gridded, 2).data

            # inverse (spectral→grid) into a fresh step-2 view each call → only ONE graph captured
            SpeedyTransforms.clear_fourier_graph_cache!()
            for _ in 1:4
                transform!(get_step(gridded, 2), specs, S)
            end
            @test sum(length(c.inverse_execs) for c in values(SpeedyTransforms.GRAPH_CACHES); init = 0) == 1

            # forward (grid→spectral) from a fresh step-2 view each call → only ONE graph captured
            SpeedyTransforms.clear_fourier_graph_cache!()
            for _ in 1:4
                transform!(specs, get_step(gridded, 2), S)
            end
            @test sum(length(c.forward_execs) for c in values(SpeedyTransforms.GRAPH_CACHES); init = 0) == 1

            SpeedyTransforms.clear_fourier_graph_cache!()
        end
    end
end
