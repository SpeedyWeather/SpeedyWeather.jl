# Shared test logic for CUDA-graphs tests, called from cuda_graphs.jl.
# HIP graphs are not implemented for AMDGPU (see SpeedyTransformsAMDGPUExt.jl and
# hip_graphs.jl), so this is CUDA-only.
# GRAPH_CACHES, MAX_GRAPHS, and clear_fourier_graph_cache! live inside the CUDA extension
# module (SpeedyTransformsCUDAExt), not in the main SpeedyTransforms package — hence
# accessed via `ext.` below, not `SpeedyTransforms.`. Keeping this code extension-local
# (rather than shared via the main package) is deliberate: it broke Enzyme's CPU-only AD
# tests when these methods existed in every Julia session regardless of which GPU backend
# was loaded.

function test_gpu_graphs(ext, prefix)
    @testset "$prefix Graphs: bounded graph cache over a GPU model run" begin
        if ext !== nothing
            spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = SpeedyWeather.GPU())
            model = PrimitiveWetModel(spectral_grid)
            simulation = initialize!(model)

            ext.clear_fourier_graph_cache!()

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
            @test s2.failed == 0                # nothing fell back to the un-captured direct loop
            @test s2.total == s1.total          # extra steps add no graphs → cache is bounded
            @test s1.total > 0                  # graphs were actually captured
            @test s1.maxlen < ext.MAX_GRAPHS    # cache stays under the cap

            ext.clear_fourier_graph_cache!()
        end
    end

    # Focused regression test for the graph-cache keying. The time stepping fetches each
    # transform's grid field via a per-step view (`get_step` → `field_view`). On GPU every
    # such view is a FRESH array wrapper aliasing the same device memory, so the cache must
    # key on the stable device pointer — not the churning wrapper identity — or it captures a
    # new graph every step and grows without bound.
    @testset "$prefix Graphs: per-step views of one buffer reuse a single graph" begin
        if ext !== nothing
            spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = SpeedyWeather.GPU())
            S = SpectralTransform(spectral_grid)
            nlayers = spectral_grid.nlayers

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

    # `add=true` (accumulate onto the field instead of overwriting; used by the Enzyme adjoint
    # rules). Two hazards, both regressions guarded here:
    #  1. A captured graph bakes in overwrite-vs-accumulate, so the two modes must NOT share a
    #     graph — the cache is therefore keyed by (field buffer, add). Sharing would silently
    #     give wrong results.
    #  2. `run_graph!` warms up and then captures; capture only RECORDS (it does not execute), so
    #     the warmup alone produces the call's result. Launching the graph as well would apply the
    #     work twice — invisible for an overwriting loop, but a double-accumulate for `add=true`.
    @testset "$prefix Graphs: add=true accumulates exactly once and uses its own graph" begin
        if ext !== nothing
            spectral_grid = SpectralGrid(; trunc = 15, nlayers = 4, architecture = SpeedyWeather.GPU())
            S = SpectralTransform(spectral_grid)
            nlayers = spectral_grid.nlayers
            specs = rand(ComplexF32, spectral_grid.spectrum, nlayers)
            field = zeros(Float32, spectral_grid.grid, nlayers)
            scratch = S.scratch_memory

            n_inverse() = sum(length(c.inverse_execs) for c in values(ext.GRAPH_CACHES); init = 0)

            ext.clear_fourier_graph_cache!()
            SpeedyTransforms._transform_grid!(field, specs, scratch, S, false)    # populate the fourier scratch

            SpeedyTransforms._fourier!(field, scratch.north, scratch.south, S)                # overwrite → reference
            once = Array(copy(field.data))
            n_overwrite = n_inverse()

            SpeedyTransforms._fourier!(field, scratch.north, scratch.south, S; add = true)    # accumulate on same buffer
            twice = Array(copy(field.data))

            @test n_inverse() > n_overwrite                 # add=true captured its OWN graph (not shared)
            @test twice ≈ 2 .* once                         # applied exactly once, not twice (no double-add)

            # replaying both cached modes stays correct and captures nothing new
            n_before = n_inverse()
            SpeedyTransforms._fourier!(field, scratch.north, scratch.south, S)
            @test Array(field.data) ≈ once
            SpeedyTransforms._fourier!(field, scratch.north, scratch.south, S; add = true)
            @test Array(field.data) ≈ 2 .* once
            @test n_inverse() == n_before

            ext.clear_fourier_graph_cache!()
        end
    end
end
