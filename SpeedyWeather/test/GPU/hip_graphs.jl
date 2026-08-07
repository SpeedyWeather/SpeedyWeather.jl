# HIP graphs are not implemented for AMDGPU (see SpeedyTransformsAMDGPUExt.jl for why);
# `gpu_graphs = true` must fall back to the generic Fourier path and warn exactly once.
ext = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsAMDGPUExt)

@testset "HIP Graphs: not implemented, falls back to the generic path with a warning" begin
    if ext !== nothing
        spectral_grid = SpectralGrid(; trunc = 31, nlayers = 8, architecture = SpeedyWeather.GPU())
        field = rand(Float32, spectral_grid.grid, spectral_grid.nlayers)

        S_off = SpectralTransform(spectral_grid; gpu_graphs = false)
        spec_off = transform(field, S_off)

        S_on = SpectralTransform(spectral_grid)     # gpu_graphs = true by default
        ext.GPU_GRAPHS_WARNED[] = false
        @test_logs (:warn, r"not implemented") transform(field, S_on)   # warns once
        spec_on = @test_logs transform(field, S_on)                     # then stays silent

        @test Array(spec_on.data) ≈ Array(spec_off.data)   # same (unaccelerated) result either way
    end
end
