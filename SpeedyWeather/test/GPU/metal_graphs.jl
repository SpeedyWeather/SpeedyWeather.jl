# Tests for the MPSGraph-fused batched Fourier transform (Metal's `gpu_graphs` path),
# analogous to cuda_graphs.jl/hip_graphs.jl + gpu_graphs_shared.jl but Metal-specific: Metal
# doesn't share the CUDA/HIP graph-capture machinery in `gpu_graphs_common.jl`
# (`GRAPH_CACHES`/`MAX_GRAPHS`/`clear_fourier_graph_cache!`) — it has its own caches
# (`FOURIER_GRAPH_CACHES`/`FOURIER_BUFFER_CACHES`/`clear_metal_fourier_cache!`) defined in
# `SpeedyTransformsMetalExt.jl`.
#
# Unlike the CUDA/AMDGPU tests, this one does NOT gate on `SpeedyTransforms.default_gpu_graphs`:
# `default_gpu_graphs(::Metal.MetalBackend)` is currently forced to `false` (see
# `SpeedyTransformsMetalExt.jl`, "Disable gpu_graphs for debugging with Metal GPU CI pipeline"),
# but whether it's safe to re-enable is exactly what's unknown — these tests exist to find out,
# so they always run `gpu_graphs = true` explicitly when the Metal extension is loaded.

function test_metal_graphs(ext)
    # Standalone grid list (not `grid_list` from spectral_transform.jl -- that file is commented
    # out of runtests.jl right now to save CI time while iterating on the synchronization fix
    # below; restore the `include` and drop this once other Metal GPU tests are re-enabled).
    grid_list = [FullGaussianGrid, OctahedralGaussianGrid, OctahedralClenshawGrid]

    @testset "Metal Graphs: fourier_batched equivalence (gpu_graphs on vs off)" begin
        if ext !== nothing
            @testset for Grid in grid_list
                spectral_grid = SpectralGrid(; truncation = 32, nlayers = 8, Grid, architecture = SpeedyWeather.GPU(), dealiasing = 3)
                field = rand(Float32, spectral_grid.grid, spectral_grid.nlayers)
                coeffs = rand(ComplexF32, spectral_grid.spectrum, spectral_grid.nlayers)

                # reference: the plain per-ring path (gpu_graphs disabled, today's Metal default)
                S_off = SpectralTransform(spectral_grid; gpu_graphs = false)
                spec_off = transform(field, S_off)       # grid -> spectral
                grid_off = transform(coeffs, S_off)       # spectral -> grid

                # MPSGraph-fused path, explicitly enabled
                ext.clear_metal_fourier_cache!()
                S_on = SpectralTransform(spectral_grid; gpu_graphs = true)
                spec_on = transform(field, S_on)
                grid_on = transform(coeffs, S_on)

                @test Array(spec_on.data) ≈ Array(spec_off.data) rtol = sqrt(eps(Float32))
                @test Array(grid_on.data) ≈ Array(grid_off.data) rtol = sqrt(eps(Float32))

                # a graph and its buffer cache were actually built for this transform
                @test haskey(ext.FOURIER_GRAPH_CACHES, S_on)
                @test haskey(ext.FOURIER_BUFFER_CACHES, S_on)

                # replaying into the same buffers across repeated calls stays correct and reuses
                # the cached graph/buffers (no growth: still exactly one entry, keyed by K=nlayers)
                n_before = length(ext.FOURIER_GRAPH_CACHES[S_on])
                spec_repeat = similar(spec_on)
                for _ in 1:3
                    transform!(spec_repeat, field, S_on)
                end
                @test Array(spec_repeat.data) ≈ Array(spec_off.data) rtol = sqrt(eps(Float32))
                @test length(ext.FOURIER_GRAPH_CACHES[S_on]) == n_before

                ext.clear_metal_fourier_cache!()
            end
        end
    end

    @testset "Metal Graphs: full PrimitiveWetModel run with gpu_graphs=true" begin
        if ext !== nothing
            spectral_grid = SpectralGrid(; truncation = 32, nlayers = 8, architecture = SpeedyWeather.GPU())
            spectral_transform = SpectralTransform(spectral_grid; gpu_graphs = true)
            model = PrimitiveWetModel(spectral_grid; spectral_transform)
            simulation = Simulation(model)
            run!(simulation, steps = 5)

            @test simulation.model.feedback.nans_detected == false

            ext.clear_metal_fourier_cache!()
        end
    end
    return nothing
end

ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsMetalExt)
test_metal_graphs(ext)
