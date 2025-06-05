import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using KernelAbstractions
using CUDA
using Printf

# Function to create test data for a given size
function setup_test_data(L, M, N, arch)
    NF = Float32
    vor = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    div = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    U = similar(vor)
    V = similar(vor)
    S = SpectralTransform(vor)
    return U, V, vor, div, S
end

# Function to run benchmark for a given size and implementation
function run_benchmark(L, M, N, arch, implementation=:original)
    U, V, vor, div, S = setup_test_data(L, M, N, arch)
    
    if implementation == :original
        return @benchmark CUDA.@sync SpeedyWeather.SpeedyTransforms.UV_from_vordiv_KA!($U, $V, $vor, $div, $S)
    elseif implementation == :split
        return @benchmark CUDA.@sync SpeedyWeather.SpeedyTransforms.UV_from_vordiv_KA_split!($U, $V, $vor, $div, $S)
    else
        return @benchmark SpeedyWeather.SpeedyTransforms.UV_from_vordiv!($U, $V, $vor, $div, $S)
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

# Function to run benchmarks with different architectures for each implementation
function run_arch_benchmarks(std_arch, ka_arch, split_arch)
    # Print header
    println("\nBenchmarking UV_from_vordiv implementations")
    println("Standard implementation architecture: ", typeof(std_arch))
    println("KA implementation architecture: ", typeof(ka_arch))
    println("KA split implementation architecture: ", typeof(split_arch))
    println("\nSize (L×M×N)      standard (", typeof(std_arch), ")    KA original (", typeof(ka_arch), ")    KA split (", typeof(split_arch), ")    Speedup (KA/std)    Speedup (split/KA)")
    println("-" ^ 120)
    
    # Run benchmarks for each size
    for (L, M, N) in sizes
        # Run all implementations on their respective architectures
        b_std = run_benchmark(L, M, N, std_arch, :standard)      # Standard UV_from_vordiv!
        b_ka = run_benchmark(L, M, N, ka_arch, :original)        # KA version UV_from_vordiv_KA!
        b_split = run_benchmark(L, M, N, split_arch, :split)     # Split kernel version UV_from_vordiv_KA_split!
        
        # Calculate median times in milliseconds
        t_std = median(b_std.times) / 1e6  # Convert ns to ms
        t_ka = median(b_ka.times) / 1e6
        t_split = median(b_split.times) / 1e6
        
        # Calculate speedups
        speedup_ka_std = t_std / t_ka
        speedup_split_ka = t_ka / t_split
        
        # Print results
        @printf("%3d×%3d×%2d     %8.3f ms    %8.3f ms    %8.3f ms      %5.2fx           %5.2fx\n",
                L, M, N, t_std, t_ka, t_split, speedup_ka_std, speedup_split_ka)
    end
end

# CPU vs CPU vs CPU benchmarks
println("\n==== CPU vs CPU vs CPU BENCHMARKS ====")
run_arch_benchmarks(SpeedyWeather.CPU(), SpeedyWeather.CPU(), SpeedyWeather.CPU())

# Run GPU benchmarks if available
if CUDA.functional()
    # CPU vs GPU vs GPU benchmarks
    println("\n==== CPU vs GPU vs GPU BENCHMARKS ====")
    run_arch_benchmarks(SpeedyWeather.CPU(), SpeedyWeather.GPU(), SpeedyWeather.GPU())
else
    println("\nCUDA GPU not available for benchmarking")
end
