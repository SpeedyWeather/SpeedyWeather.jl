# Test if the model runs on Reactant with GPU as well

using Reactant, CUDA # we need CUDA even on other devices for the kernel raising

Reactant.set_default_backend("gpu")

arch = SpeedyWeather.ReactantDevice()
spectral_grid = SpectralGrid(architecture = arch, nlayers = 1)
spectral_transform = MatrixSpectralTransform(spectral_grid)
model = BarotropicModel(
    spectral_grid = spectral_grid,
    spectral_transform = spectral_transform,
    output = nothing,
    feedback = nothing,
)
simulation = initialize!(model)
run!(simulation, steps = 10)

@test !any(isnan, on_architecture(SpeedyWeather.CPU(), simulation.variables.prognostic.vorticity))

# also test the PrimitiveWetModel
spectral_grid = SpectralGrid(architecture = arch)
spectral_transform = MatrixSpectralTransform(spectral_grid)
model = PrimitiveWetModel(spectral_grid; spectral_transform, feedback = nothing, output = nothing)
simulation = initialize!(model)
run!(simulation, steps = 10)

@test !any(isnan, on_architecture(SpeedyWeather.CPU(), simulation.variables.prognostic.vorticity))
