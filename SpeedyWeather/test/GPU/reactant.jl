# Reactant on the GPU: the BarotropicModel steps without NaNs. Not a correctness test (too expensive
# here, see test/reactant for the correctness tests on Reactant's CPU backend). A standard model run,
# Reactant raises the kernels by default (raise = true) in `initialize!` (`@maybe_jit`) and `run!`

using Reactant, CUDA # we need CUDA even on other devices for the kernel raising

Reactant.set_default_backend("gpu")

@testset "Reactant GPU: BarotropicModel steps without NaNs" begin
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice(), nlayers = 1)
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = BarotropicModel(spectral_grid; spectral_transform, feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 10)

    @test !any(isnan, on_architecture(SpeedyWeather.CPU(), simulation.variables.prognostic.vorticity))
end

#= # TODO: Reactative test, the compiled_call just takes too long in the CI with consumer GPUs :(
# with both spectral transforms (the FFT-based one uses cuFFT through XLA on the GPU)
@testset "Reactant GPU: PrimitiveWetModel with $Transform steps without NaNs" for Transform in (MatrixSpectralTransform, ) # SpectralTransform
    spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice())
    spectral_transform = Transform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform, feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 10)

    (; prognostic) = simulation.variables
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
        @test !any(isnan, on_architecture(SpeedyWeather.CPU(), getproperty(prognostic, name)))
    end
end
=# 