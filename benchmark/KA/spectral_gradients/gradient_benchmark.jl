# Benchmark the performance of ∇²! and ∇²_KA!, comparing the old non-kernel version with the new KA version

import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using CUDA
using KernelAbstractions
using Printf

const DEFAULT_RADIUS = 1

# old version non-KA
function ∇_old!(
    dpdx::LowerTriangularArray,     # Output: zonal gradient
    dpdy::LowerTriangularArray,     # Output: meridional gradient
    p::LowerTriangularArray,        # Input: spectral coefficients
    S::SpectralTransform;           # includes precomputed arrays
    radius = DEFAULT_RADIUS,        # scale with radius if provided, otherwise unit sphere
)
    (; grad_y1, grad_y2) = S
    @boundscheck SpeedyWeather.SpeedyTransforms.ismatching(S, p) || throw(DimensionMismatch(S, p))

    # maximum degree l, order m of spherical harmonics (1-based)
    lmax, mmax = size(p, OneBased, as=Matrix)

    for k in eachmatrix(dpdx, dpdy, p)      # also performs size checks
        lm = 0
        @inbounds for m in 1:mmax-1         # 1-based l, m, skip last column

            # DIAGONAL (separated to avoid access to l-1, m which is above the diagonal)
            lm += 1

            dpdx[lm, k] = (m-1)*im*p[lm, k]         # zonal gradient: d/dlon = *i*m
            dpdy[lm, k] = grad_y2[lm]*p[lm+1, k]    # meridional gradient: p[lm-1]=0 on diagonal
        
            # BELOW DIAGONAL (all terms)
            for l in m+1:lmax-1                     # skip last row
                lm += 1
                dpdx[lm, k] = (m-1)*im*p[lm, k]
                dpdy[lm, k] = grad_y1[lm]*p[lm-1, k] + grad_y2[lm]*p[lm+1, k]
            end

            # LAST ROW (separated to avoid out-of-bounds access to lmax+1
            lm += 1
            dpdx[lm, k] = (m-1)*im*p[lm, k]
            dpdy[lm, k] = grad_y1[lm]*p[lm-1, k]    # only first term from 2nd last row
        end

        # LAST COLUMN
        @inbounds begin
            lm += 1                                 # second last row
            dpdx[lm, k] = (mmax-1)*im*p[lm, k]
            dpdy[lm, k] = grad_y2[lm]*p[lm+1, k]    # only 2nd term

            lm += 1                                 # last row
            dpdx[lm, k] = (mmax-1)*im*p[lm, k]
            dpdy[lm, k] = grad_y1[lm]*p[lm-1, k]    # only 1st term
        end
    end

    # 1/radius factor if not unit sphere
    if radius != 1
        R⁻¹ = inv(radius)
        dpdx .*= R⁻¹
        dpdy .*= R⁻¹
    end

    return dpdx, dpdy
end


# Function to create test data for a given size
function setup_test_data(L, M, N, arch)
    NF = Float32
    alms = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
    alms_dx_out = similar(alms)
    alms_dy_out = similar(alms)
    S = SpectralTransform(alms)
    return alms_dx_out, alms_dy_out, alms, S
end

# Function to run benchmark for a given size and implementation
function run_benchmark(L, M, N, arch, use_ka=true)
    alms_dx_out, alms_dy_out, alms, S = setup_test_data(L, M, N, arch)
    
    if use_ka
        return @benchmark CUDA.@sync SpeedyWeather.SpeedyTransforms.∇!($alms_dx_out, $alms_dy_out, $alms, $S)
    else
        return @benchmark ∇_old!($alms_dx_out, $alms_dy_out, $alms, $S)
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
    println("\nBenchmarking ∇ implementations")
    println("Standard implementation architecture: ", typeof(std_arch))
    println("KA implementation architecture: ", typeof(ka_arch))
    println("\nSize (L×M×N)      old ∇! (", typeof(std_arch), ")    new (KA) ∇! (", typeof(ka_arch), ")     Speedup")
    println("-" ^ 100)
    
    # Run benchmarks for each size
    for (L, M, N) in sizes
        # Run both implementations on their respective architectures
        b_std = run_benchmark(L, M, N, std_arch, false)  # Standard ∇!
        b_ka = run_benchmark(L, M, N, ka_arch, true)     # KA version ∇_KA!
        
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

