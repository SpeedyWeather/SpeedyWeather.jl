# Test if the model runs on Reactant with GPU as well 

using Reactant, CUDA # we need CUDA even on other devices for the kernel raising 

Reactant.set_default_backend("gpu")

arch = SpeedyWeather.ReactantDevice()
spectral_grid = SpectralGrid(architecture=arch, nlayers=1)
spectral_transform = MatrixSpectralTransform(spectral_grid)
model = BarotropicModel(spectral_grid=spectral_grid, spectral_transform=spectral_transform)
simulation = initialize!(model)
run!(simulation, nsteps=10)

# TODO: fix GPU reactant in follow up when initial conditions are easier to handle
@test_broken !any(isnan.(simulation.prognostic_variables.vor))


