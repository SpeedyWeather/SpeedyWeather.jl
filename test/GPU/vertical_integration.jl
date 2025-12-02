# The vertical_integration! and vertical_velocity! function are
# some of the few functions that have distinct GPU version. 
# Here, we just test that they are identical to the CPU version
# The CPU version is in turn tested more properly in the regular unit test

@testset "Vertical Integration" begin
    # TODO: Add actual tests for vertical integration

    spectral_grid_cpu = SpectralGrid(architecture=CPU())
    spectral_grid_gpu = SpectralGrid(architecture=GPU())

    model_cpu = PrimitiveWetModel(spectral_grid_cpu)
    model_gpu = PrimitiveWetModel(spectral_grid_gpu)
end