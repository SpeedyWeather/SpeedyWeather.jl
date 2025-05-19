# Benchmark the performance of ∇²! and ∇²_KA!, comparing the old non-kernel version with the new KA version

import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using CUDA
using KernelAbstractions
using Printf

# this is the old version that's not there anymore 
function old_∇²!(
    ∇²alms::LowerTriangularArray,   # Output: (inverse) Laplacian of alms
    alms::LowerTriangularArray,     # Input: spectral coefficients
    S::SpectralTransform;           # precomputed eigenvalues
    add::Bool=false,                # add to output array or overwrite
    flipsign::Bool=false,           # -∇² or ∇²
    inverse::Bool=false,            # ∇⁻² or ∇²
    radius = DEFAULT_RADIUS,        # scale with radius if provided, otherwise unit sphere
)
    @boundscheck ismatching(S, ∇²alms) || throw(DimensionMismatch(S, ∇²alms))

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
]

# Choose architecture based on CUDA availability
arch = CUDA.functional() ? SpeedyWeather.CUDAGPU() : SpeedyWeather.CPU()

# Print header
println("\nBenchmarking ∇² implementations")
println("Architecture: ", typeof(arch))
println("\nSize (L×M×N)      old ∇²! (CPU)          new (KA) ∇²! (", typeof(arch), ")     Speedup")
println("-" ^ 70)

# Run benchmarks for each size
for (L, M, N) in sizes
    # Run both implementations
    b_std = run_benchmark(L, M, N, arch, false)  # Standard ∇²!
    b_ka = run_benchmark(L, M, N, arch, true)    # KA version ∇²_KA!
    
    # Calculate median times in milliseconds
    t_std = median(b_std.times) / 1e6  # Convert ns to ms
    t_ka = median(b_ka.times) / 1e6
    speedup = t_std / t_ka
    
    # Print results
    @printf("%3d×%3d×%2d     %8.3f ms      %8.3f ms         %5.2fx\n",
            L, M, N, t_std, t_ka, speedup)
end
println("Running benchmarks on ", typeof(arch))

