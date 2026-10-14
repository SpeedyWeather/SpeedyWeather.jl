# Reactant on the GPU: the PrimitiveWetModel steps without NaNs. Not a correctness test (too expensive
# here, see test/reactant for the correctness tests on Reactant's CPU backend). A standard model run,
# Reactant raises the kernels by default (raise = true) in `initialize!` (`@maybe_jit`) and `run!`

using Reactant, CUDA # we need CUDA even on other devices for the kernel raising

Reactant.set_default_backend("gpu")

@testset "Reactant GPU: PrimitiveWetModel steps without NaNs" begin
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice())
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform, feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 10)

    (; prognostic) = simulation.variables
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
        @test !any(isnan, on_architecture(SpeedyWeather.CPU(), getproperty(prognostic, name)))
    end
end
