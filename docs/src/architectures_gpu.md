# GPU & Architectures

!!! warning "Work in progress"
    The GPU support of SpeedyWeather.jl is still work in progress and some parts of this documentation might not be always updated to the latest state. We will extend this documentation over time. Don't hesitate to contact us via GitHub issues or mail when you have questions or want to collaborate.

Most of SpeedyWeather.jl supports GPU acceleration. All of our models can run GPUs, however as our development of this is still very recent, there still might be issues with the GPU models. If you encounter any, please file a GitHub issue. Our development focuses on CUDA GPUs, but other architectures are thinkable in the future as well, as our approach relies on the device agnostic `KernelAbstractions.jl`. An experimental port to AMD GPUs using the `AMDGPU` package is available but AMD-specific performance optimizations are not implemented yet. The SpeedyWeather.jl submodule `Architectures` encodes all the information of the device we run our models on. In order to initialize a model on a GPU, we need to load the `CUDA` or `AMDGPU` package and pass the architecture to the model constructor. For example, to initialize a barotropic model on a GPU, we can do the following:  

```julia
using SpeedyWeather, CUDA # For AMD GPUs, replace `CUDA` with `AMDGPU`
architecture = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(trunc=41, nlayers=1, architecture=architecture)           

model = PrimitiveWetModel(spectral_grid=spectral_grid)
simulation = initialize!(model)
run!(simulation, period=Day(10))
```

## Architectures Utilities 

In order to easily transfer our structures between CPU (e.g. for plotting and output) and GPU, we have the following utilities that can make use of the `architecture` object defined above and the `on_architecture` function, e.g. as follows: 

```julia
using SpeedyWeather, CUDA # For AMD GPUs, replace `CUDA` with `AMDGPU`
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

Be aware that directly calling e.g. `CuArray`, `ROCArray` or `adapt` on the data structres is not recommended, as it can lead to unexpected behavior, e.g. mismatching internal architecture representations when launching kernels and other operations. Please use the `on_architecture` function instead for all transfer between devices. 

## Multithreading 

Our implementation of the model using KernelAbstractions.jl, also enables an easy multithreading of the model on CPU. As soon as you start Julia with more than one thread (e.g. `julia --threads 4`) SpeedyWeather will make use of all threads available. 

## Benchmarks 

More to follow...
