# GPU and Architectures

!!! warning "Work in progress"
    The GPU support of SpeedyWeather.jl is still work in progress and some parts of this documentation
    might not be always updated to the latest state. We will extend this documentation over time.
    Don't hesitate to contact us via GitHub issues or mail when you have questions or want to collaborate.

Most of SpeedyWeather.jl supports GPU acceleration. All of our models can run GPUs, however as our
development of this is still very recent, there still might be issues with the GPU models and do
not consider them fully optimized yet. If you encounter any issues, please report them via GitHub issues.
Our development focuses on CUDA GPUs, but other architectures are thinkable in the future as well,
as our approach relies on the device agnostic `KernelAbstractions.jl`. An experimental port to AMD
GPUs using the `AMDGPU` package is available; some AMD-specific performance optimizations exist
(e.g. GPU graphs, see below) but are not enabled by default yet, pending broader hardware verification.
There is also support for Apple's `Metal` package, letting you run SpeedyWeather.jl on the GPU of
Apple Silicon Macs; this is convenient for local development and small experiments, but Metal has
no double-precision (Float64) support, so use `NF = Float32` (the default) throughout when using it.
The SpeedyWeather.jl submodule `Architectures` encodes all the information of the device we run our models on.
In order to initialize a model on a GPU, we need to load the `CUDA`, `AMDGPU`, or `Metal` package and
pass the architecture to the model constructor. For example, to initialize a barotropic model on a
GPU, we can do the following:  

```julia
using SpeedyWeather, CUDA # For AMD GPUs use `AMDGPU`, for Apple Silicon use `Metal`
architecture = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(truncation = 64, nlayers = 8, architecture = architecture)           

model = PrimitiveWetModel(spectral_grid=spectral_grid)
simulation = initialize!(model)
run!(simulation, period=Day(10))
```

## GPU Graphs

On GPU, the batched Fourier transform (part of the spectral transform) can be accelerated with
GPU graphs (CUDA graphs / HIP graphs). With reduced ring grids used, rings / latitudes have different
numbers of grid points per latitude. 
This means that the Fourier transform can't be computed as just one large batched FFT. 
Launching many FFTs for a single transform causes significant launch overheads on GPU.
The GPU graphs eliminate that by capturing and replaying these launches. 
The fused gather/FFT/scatter loop is captured once and
replayed with a single launch on every subsequent call, instead of re-issuing many small
kernel/FFT launches per step. This mainly helps at lower resolutions, where the transform is
launch-bound rather than compute-bound.

It is controlled by the `gpu_graphs` keyword of `SpectralTransform`, whose default is
backend-dependent (`SpeedyTransforms.default_gpu_graphs`) rather than always on. Models don't
expose `gpu_graphs` directly, so enable or disable it by constructing your own transform and
passing it to the model:

```julia
spectral_transform = SpectralTransform(spectral_grid; gpu_graphs = true)
model = PrimitiveWetModel(spectral_grid; spectral_transform)
```

On **CUDA**, GPU graphs are enabled by default.

On **AMDGPU**, the HIP-graphs implementation mirrors the CUDA one but is **disabled by default**.
It has been verified stable on datacenter/CDNA hardware (e.g. LUMI), but caused failures on the
consumer RDNA cards our CI currently runs on, so it isn't trusted as a default yet. To try it
explicitly on hardware you trust, pass `gpu_graphs = true` as shown above.

On **Metal**, `gpu_graphs` is **disabled by default** and not currently supported; the Fourier
transform instead batches north/south ring pairs into single Metal Performance Shaders Graph
(`MPSGraph`) calls to reduce CPU/GPU communication overhead.

## Architectures Utilities 

In order to easily transfer our structures between CPU (e.g. for plotting and output) and GPU,
we have the following utilities that can make use of the `architecture` object defined above
and the `on_architecture` function, e.g. as follows: 

```julia
using SpeedyWeather, CUDA # For AMD GPUs use `AMDGPU`, for Apple Silicon use `Metal`
nlat_half = 6
arch_cpu = SpeedyWeather.CPU()
arch_gpu = SpeedyWeather.GPU()

grid_cpu = HEALPixGrid(nlat_half, arch_cpu)
grid_gpu = on_architecture(arch_gpu, grid_cpu)

field_cpu = rand(grid_cpu)
field_gpu = on_architecture(arch_gpu, field_cpu)

spectrum_cpu = Spectrum(truncation=43, architecture=arch_cpu)
spectrum_gpu = on_architecture(arch_gpu, spectrum_cpu)

spec_cpu = rand(spectrum_cpu)
spec_gpu = on_architecture(arch_gpu, spec_cpu)
```

Be aware that directly calling e.g. `CuArray`, `ROCArray`, `MtlArray` or `adapt` on the data structres
is not recommended, as it can lead to unexpected behavior, e.g. mismatching internal architecture
representations when launching kernels and other operations. Please use the `on_architecture`
function instead for all transfer between devices. See also the [GPU](@ref) section of
`LowerTriangularArrays` for details on the GPU support of spectral arrays, and 
[SpeedyTransforms](@ref) for spectral transforms on GPU.

## Multithreading 

Our implementation of the model using KernelAbstractions.jl, also enables an easy multithreading
of the model on CPU. As soon as you start Julia with more than one thread (e.g. `julia --threads 4`)
SpeedyWeather will make use of all threads available. 

## Reactant Support 

!!! warning "Work in progress"
    The Reactant support of SpeedyWeather.jl is still work in progress, incomplete and considered experimental.

SpeedyWeather.jl is compatible with the MLIR/XLA optimizer Reactant.jl that enables to run an optimized OpenXLA accelearted version of our model (this is the same library that powers e.g. JAX numerics). It also enables to run the model on CPU, different GPUs and even TPU architectures. 
All models (`BarotropicModel`, `ShallowWaterModel`, `PrimitiveDryModel`, `PrimitiveWetModel`)
run with Reactant, use the `ReactantDevice` architecture and set your device within Reactant with
`Reactant.set_default_backend`:

```julia
using Reactant, CUDA   # CUDA is needed for Reactant to compile KernelAbstractions kernels
using SpeedyWeather

spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice(), nlayers = 8)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period = Day(10))
```

With `ReactantDevice` the models use the `SpectralTransform` by default. On Reactant it doesn't use
FFT plans or scratch memory; both of its stages are matrix multiplies that XLA compiles together with
the rest of the time step:

- The Fourier transform is a truncated discrete Fourier transform (`fourier_dft = true`, the default up
  to truncation 256, `SpeedyTransforms.DFT_MAX_TRUNCATION`). It computes only the orders that the
  Legendre transform uses, as one matrix multiply per direction batched over all latitude rings. Its
  matrices need memory that grows with the cube of the truncation: 1.4 MB at T32, 617 MB at T256 in
  Float32. With `fourier_dft = false` (the default above truncation 256) it uses one FFT per group of
  latitude rings of equal length instead. XLA can't put FFTs into its CUDA graphs, though, which splits
  a compiled time step on GPUs into many separately launched kernels.
- The Legendre transform is one matrix multiply per direction, batched over the spherical harmonic
  order.

`run!` compiles the time loop with Reactant into programs that run many time steps at once
(a traced loop, see `SpeedyWeather.time_steps!`), as calling a compiled function from Julia costs
more time than a time step takes on the GPU at lower resolutions. Output can't be written from
inside a compiled program, so the time steps between two outputs are one program and the output
is written between the programs. Without output (`run!(simulation, output = false)`, the default)
the whole run is one program. With output, the spin-up steps of the time stepper at the start
and a remainder (if the run isn't a whole number of output intervals) use a compiled single time
step, so at most two programs are compiled per `run!`. A program of several time steps takes about
twice as long to compile as a single time step. The progress meter and the check for NaNs (with
`Feedback(debug = true)`, the default) are updated after every program, not every time step. SpeedyWeather compiles with `raise = true` (and without the optimization patterns listed in `SpeedyWeather.REACTANT_EXCLUDED_PASSES`, workarounds for bugs in Reactant; compile your own programs with `SpeedyWeather.reactant_compile_options()` to get the same options): Reactant
raises the KernelAbstractions kernels to StableHLO, the operations of the XLA compiler, which on
GPUs is needed for correct results. Not all model components are compatible yet
(e.g. particle advection). The setups we use in the CI to compare the results with Reactant to
those on the CPU are in `SpeedyWeather/test/reactant/`, on CUDA GPUs our GPU tests check that a
`BarotropicModel` runs without NaNs with both transforms (`SpeedyWeather/test/GPU/reactant.jl`).

### Precision of matrix products on GPUs

On NVIDIA GPUs with TF32 support (Ampere and newer, e.g. A40, A100, H100) XLA runs Float32 matrix
products in TF32 by default, which rounds the inputs to a 10-bit mantissa: relative errors of about
1e-4 per product instead of about 1e-7 in Float32. This concerns the matrix products of the
`MatrixSpectralTransform`, the Legendre transform of the `SpectralTransform` (and its Fourier
transform as truncated DFT, see above) and those that XLA
creates from loops in raised kernels (e.g. the vertical
matrix-vector products of the semi-implicit time stepping).

The model runs stably with TF32 but the results deviate considerably from those in Float32 (e.g.
on the CPU or with CUDA.jl). For a `PrimitiveWetModel` (T31, 8 layers, `MatrixSpectralTransform`) on
an NVIDIA A40 the maximum relative differences to the CPU after 10 time steps are

| Variable         | TF32 (default) | `HIGHEST` (see below) |
| ---------------- | -------------- | --------------------- |
| vorticity        | 1.4e-2         | 8.0e-5                |
| divergence       | 3.9e-2         | 2.6e-4                |
| temperature      | 5.0e-4         | 1.8e-6                |
| humidity         | 8.4e-3         | 1.9e-5                |
| surface pressure | 4.7e-6         | 2.7e-8                |

and after the first time step the divergence even differs by 48% with TF32. So use Float32
precision whenever the results should agree with Float32 results, e.g. for comparisons, tests or
gradients (see [Differentiating with Reactant](@ref)). Whether TF32 is acceptable for your
application, e.g. for long climate simulations, we haven't evaluated.

For Float32 precision wrap everything that compiles, i.e. `initialize!` and `run!`, in
`Reactant.with_config(dot_general_precision = Reactant.PrecisionConfig.HIGHEST)`:

```julia
using Reactant, CUDA
using SpeedyWeather

spectral_grid = SpectralGrid(architecture = SpeedyWeather.ReactantDevice(), nlayers = 8)
model = PrimitiveWetModel(spectral_grid)

simulation = Reactant.with_config(dot_general_precision = Reactant.PrecisionConfig.HIGHEST) do
    simulation = initialize!(model)
    run!(simulation, period = Day(10))
    simulation
end
```

With the default `SpectralTransform` this makes a `PrimitiveWetModel` time step (8 layers) on an H100
up to 16% slower up to T128 and 35% slower at T256 (with the `MatrixSpectralTransform` 6% at T32 to
2.3x at T128, which gave the `HIGHEST` column in the table above). In the transform alone, Float32 instead
of TF32 reduces the error from 2e-4 to 4e-4 to about 1e-7 to 9e-7, the level of the Float32 transform on
the CPU. It applies to the matrix products that Reactant creates itself, e.g. those
of the spectral transform, not (yet) to those that XLA creates from loops in raised kernels. To have
those in Float32 too, also set the environment variables `XLA_FLAGS=--xla_gpu_enable_triton_gemm=false`
(XLA's Triton GEMMs ignore the next setting) and `NVIDIA_TF32_OVERRIDE=0` (for cuBLAS) before starting
Julia. In our tests this didn't change the results of the model noticeably, but it made the time step
1.6x (T63) to 2.3x (T31) slower than with TF32.

## Benchmarks

See the [Benchmarks](@ref benchmarks_page) page for SYPD numbers across architectures
(CPU, GPU, Reactant), one section per architecture, and comparison figures of
the PrimitiveWet resolution sweep.
