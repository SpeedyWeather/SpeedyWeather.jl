# GPU & Architectures

!!! warning "Work in progress"
    The GPU support of SpeedyWeather.jl is still work in progress and some parts of this documentation might not be always updated to the latest state. We will extend this documentation over time. Don't hesitate to contact us via GitHub issues or mail when you have questions or want to collaborate.

Some of SpeedyWeather.jl already supports GPU acceleration, e.g. the barotropic model. Our development focuses on CUDA GPUs, but other architectures are thinkable in the future as well, as our approach relies on the device agnostic `KernelAbstractions.jl`. The SpeedyWeather.jl submodule `Architectures` encodes all the information of the device we run our models on. In order to initialize a model on a GPU, we need to load the `CUDA` package and pass the architecture to the model constructor. For example, to initialize a barotropic model on a GPU, we can do the following:  

```julia
using SpeedyWeather, CUDA 
architecture = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(trunc=41, nlayers=1, architecture=architecture)           

model = BarotropicModel(spectral_grid=spectral_grid)
CUDA.@allowscalar simulation = initialize!(model)
run!(simulation, period=Day(10))
```

Note that we need to use `CUDA.@allowscalar` here during initialization. Currently we do not yet support a fully GPU-accelerated model construction and initialization.

## Architectures Utilities 

In order to easily transfer our structures between CPU (e.g. for plotting and output) and GPU, we have the following utilities that make can make use of the `architecture` object defined above and the `on_architecture` function, e.g. as follows: 

```julia
using SpeedyWeather, CUDA 
nlat_half = 6
arch_cpu = SpeedyWeather.CPU()
arch_gpu = SpeedyWeather.GPU()

grid_cpu = HEALPixGrid(nlat_half, arch_cpu)
grid_gpu = on_architecture(arch_gpu, grid_cpu)

field_cpu = rand(grid_cpu)
field_gpu = on_architecture(arch_gpu, field_cpu)

spectrum_cpu = Spectrum(trunc=41, architecture=arch_cpu)
spectrum_gpu = on_architecture(arch_gpu, spectrum_cpu)

spec_cpu = rand(spectrum_cpu)
spec_gpu = on_architecture(arch_gpu, spec_cpu)
```

Be aware that directly calling e.g. `CuArray` or `adapt` on the data structres is not recommended, as it can lead to unexpected behavior, e.g. mismatching internal architecture representations when launching kernels and other operations. Please use the `on_architecture` function instead for all transfer between devices. 

## Benchmarks 

More to follow...
