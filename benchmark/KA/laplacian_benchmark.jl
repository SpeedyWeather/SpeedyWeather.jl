# Benchmark the performance of ∇²! and ∇²_KA!, comparing the old non-kernel version with the new KA version

import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using KernelAbstractions
using CUDA
using Printf

# this is the old version that's not there anymore 
function old_∇²!(
    ∇²alms::LowerTriangularArray,   # Output: (inverse) Laplacian of alms
    alms::LowerTriangularArray,     # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    add::Bool=false,                # add to output array or overwrite
    flipsign::Bool=false,           # -∇² or ∇²
    inverse::Bool=false,            # ∇⁻² or ∇²
    radius = 1,        # scale with radius if provided, otherwise unit sphere
)
    @boundscheck SpeedyWeather.SpeedyTransforms.ismatching(S, ∇²alms) || throw(DimensionMismatch(S, ∇²alms))

    # use eigenvalues⁻¹/eigenvalues for ∇⁻²/∇² based but name both eigenvalues
    eigenvalues = inverse ? S.eigenvalues⁻¹ : S.eigenvalues

    kernel = flipsign ? (add ? (o,a) -> (o-a) : (o, a) -> -a) : 
                        (add ? (o,a) -> (o+a) : (o, a) -> a)
    
    # maximum degree l, order m of spherical harmonics (1-based)
    lmax, mmax = size(alms, OneBased, as=Matrix)

    for k in eachmatrix(∇²alms, alms)
        lm = 0
        @inbounds for m in 1:mmax
            for l in m:lmax
                lm += 1
                ∇²alms[lm, k] = kernel(∇²alms[lm, k], alms[lm, k]*eigenvalues[l])
            end
        end
    end

    # /radius² or *radius² scaling if not unit sphere
    if radius != 1
        R_plusminus_squared = inverse ? radius^2 : inv(radius^2)
        ∇²alms .*= R_plusminus_squared
    end

    return ∇²alms
end

# Function to create test data for a given size
function setup_test_data(L, M, N, arch)
    NF = Float32
    alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    alms_out = similar(alms)
    S = SpectralTransform(alms)
    return alms_out, alms, S
end

# Function to run benchmark for a given size and implementation
function run_benchmark(L, M, N, arch, use_ka=true)
    alms_out, alms, S = setup_test_data(L, M, N, arch)
    
    if use_ka
        return @benchmark SpeedyWeather.SpeedyTransforms.∇²!($alms_out, $alms, $S)
    else
        return @benchmark old_∇²!($alms_out, $alms, $S)
    end
end

# Test sizes (L, M, N)
sizes = [
    (33, 32, 1),
    (33, 32, 8),
    (65, 64, 1),
    (65, 64, 8),
    (129, 128, 1),
    (129, 128, 8),
    (257, 256, 8),
    (513, 512, 8),
    (513, 512, 16),
    (513, 512, 32)
]

# Function to run benchmarks with two different architectures
function run_arch_benchmarks(std_arch, ka_arch)
    # Print header
    println("\nBenchmarking ∇² implementations")
    println("Standard implementation architecture: ", typeof(std_arch))
    println("KA implementation architecture: ", typeof(ka_arch))
    println("\nSize (L×M×N)      old ∇²! (", typeof(std_arch), ")    new (KA) ∇²! (", typeof(ka_arch), ")     Speedup")
    println("-" ^ 100)
    
    # Run benchmarks for each size
    for (L, M, N) in sizes
        # Run both implementations on their respective architectures
        b_std = run_benchmark(L, M, N, std_arch, false)  # Standard ∇²!
        b_ka = run_benchmark(L, M, N, ka_arch, true)     # KA version ∇²_KA!
        
        # Calculate median times in milliseconds
        t_std = median(b_std.times) / 1e6  # Convert ns to ms
        t_ka = median(b_ka.times) / 1e6
        speedup = t_std / t_ka
        
        # Print results
        @printf("%3d×%3d×%2d     %8.3f ms      %8.3f ms         %5.2fx\n",
                L, M, N, t_std, t_ka, speedup)
    end
end

# CPU vs CPU benchmarks
println("\n==== CPU vs CPU BENCHMARKS ====")
run_arch_benchmarks(SpeedyWeather.CPU(), SpeedyWeather.CPU())

# Run GPU benchmarks if available
if CUDA.functional()
    # CPU vs GPU benchmarks
    println("\n==== CPU vs GPU BENCHMARKS ====")
    run_arch_benchmarks(SpeedyWeather.CPU(), SpeedyWeather.GPU())
else
    println("\nCUDA GPU not available for benchmarking")
end

