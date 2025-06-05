# Benchmark the performance of divergence! and divergence_KA!, comparing the standard version with the KA version

import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using CUDA
using KernelAbstractions
using Printf

# Function to create test data for a given size
function setup_test_data(L, M, N, arch)
    NF = Float32
    u = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    v = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    div_out = similar(u)
    S = SpectralTransform(u)
    return div_out, u, v, S
end

# Function to run benchmark for a given size and implementation
function run_benchmark(L, M, N, arch, use_ka=true)
    div_out, u, v, S = setup_test_data(L, M, N, arch)
    
    if use_ka
        return @benchmark SpeedyWeather.SpeedyTransforms.divergence_KA!($div_out, $u, $v, $S)
    else
        return @benchmark SpeedyWeather.SpeedyTransforms.divergence!($div_out, $u, $v, $S)
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
    println("\nBenchmarking divergence implementations")
    println("Standard implementation architecture: ", typeof(std_arch))
    println("KA implementation architecture: ", typeof(ka_arch))
    println("\nSize (L×M×N)      standard divergence! (", typeof(std_arch), ")    KA divergence_KA! (", typeof(ka_arch), ")     Speedup")
    println("-" ^ 100)
    
    # Run benchmarks for each size
    for (L, M, N) in sizes
        # Run both implementations on their respective architectures
        b_std = run_benchmark(L, M, N, std_arch, false)  # Standard divergence!
        b_ka = run_benchmark(L, M, N, ka_arch, true)     # KA version divergence_KA!
        
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