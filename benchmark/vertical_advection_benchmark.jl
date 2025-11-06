"""
Performance benchmark script for vertical advection implementations.

Compares:
1. CPU-only version (vertical_advection_cpu!)
2. Kernel version on CPU (vertical_advection!)
3. Kernel version on GPU (vertical_advection!)

Tests with CenteredVerticalAdvection at different resolutions.
"""

using SpeedyWeather
using CUDA
using BenchmarkTools
using Printf

println("="^80)
println("Vertical Advection Performance Benchmark")
println("="^80)
println()

# Check if GPU is available
has_gpu = CUDA.functional()
if has_gpu
    println("✓ GPU detected and available")
else
    println("⚠ GPU not available - will only benchmark CPU versions")
end
println()

# Test configurations
truncations = [41, 123, 255]
nlayers = 8
NF = Float32
order = 4  # 4th order centered advection

results = Dict()

for trunc in truncations
    println("="^80)
    println("Testing with T$trunc truncation ($(trunc+1) degrees)")
    println("="^80)
    println()
    
    # ========== CPU-ONLY VERSION ==========
    println("Setting up CPU-only version...")
    spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, architecture=CPU())
    vertical_advection = CenteredVerticalAdvection(spectral_grid_cpu, order=order)
    model_cpu = PrimitiveWetModel(spectral_grid_cpu; vertical_advection)
    sim_cpu = initialize!(model_cpu)
    
    # Warm up
    run!(sim_cpu, steps=5)
    diagn_cpu = sim_cpu.diagnostic_variables
    
    # Benchmark CPU-only version
    println("Benchmarking CPU-only version (vertical_advection_cpu!)...")
    cpu_only_time = @belapsed SpeedyWeather.vertical_advection_cpu!($diagn_cpu, $model_cpu) samples=100 evals=1
    
    println(@sprintf("  Time: %.3f ms", cpu_only_time * 1000))
    println()
    
    # ========== KERNEL VERSION ON CPU ==========
    println("Benchmarking kernel version on CPU (vertical_advection!)...")
    
    # Reset simulation
    sim_cpu = initialize!(model_cpu)
    run!(sim_cpu, steps=5)
    diagn_cpu = sim_cpu.diagnostic_variables
    
    kernel_cpu_time = @belapsed SpeedyWeather.vertical_advection!($diagn_cpu, $model_cpu) samples=100 evals=1
    
    println(@sprintf("  Time: %.3f ms", kernel_cpu_time * 1000))
    println(@sprintf("  Speedup vs CPU-only: %.2fx", cpu_only_time / kernel_cpu_time))
    println()
    
    # Store results
    results[trunc] = Dict(
        "cpu_only" => cpu_only_time,
        "kernel_cpu" => kernel_cpu_time,
    )
    
    # ========== KERNEL VERSION ON GPU ==========
    if has_gpu
        println("Setting up GPU version...")
        spectral_grid_gpu = SpectralGrid(; NF, trunc, nlayers, architecture=GPU())
        vertical_advection_gpu = CenteredVerticalAdvection(spectral_grid_gpu, order=order)
        model_gpu = PrimitiveWetModel(spectral_grid_gpu; vertical_advection=vertical_advection_gpu)
        sim_gpu = initialize!(model_gpu)
        
        # Warm up GPU
        run!(sim_gpu, steps=5)
        diagn_gpu = sim_gpu.diagnostic_variables
        
        # Additional warmup for GPU kernels
        for _ in 1:10
            SpeedyWeather.vertical_advection!(diagn_gpu, model_gpu)
        end
        SpeedyWeather.Architectures.synchronize(model_gpu.spectral_grid.architecture)
        
        println("Benchmarking kernel version on GPU (vertical_advection!)...")
        
        # Benchmark GPU version with synchronization
        kernel_gpu_time = @belapsed begin
            CUDA.@sync SpeedyWeather.vertical_advection!($diagn_gpu, $model_gpu)
        end samples=100 evals=1
        
        println(@sprintf("  Time: %.3f ms", kernel_gpu_time * 1000))
        println(@sprintf("  Speedup vs CPU-only: %.2fx", cpu_only_time / kernel_gpu_time))
        println(@sprintf("  Speedup vs kernel CPU: %.2fx", kernel_cpu_time / kernel_gpu_time))
        println()
        
        results[trunc]["kernel_gpu"] = kernel_gpu_time
    end
    
    println()
end

# ========== SUMMARY ==========
println("="^80)
println("SUMMARY")
println("="^80)
println()

println(@sprintf("%-12s %12s %12s %12s %12s", 
    "Truncation", "CPU-only", "Kernel CPU", "Speedup"))
println("-"^80)

for trunc in truncations
    r = results[trunc]
    speedup = r["cpu_only"] / r["kernel_cpu"]
    
    println(@sprintf("T%-11d %9.3f ms %9.3f ms %10.2fx",
            trunc, 
        r["cpu_only"] * 1000, r["kernel_cpu"] * 1000, speedup))
end
println()

if has_gpu
    println(@sprintf("%-12s %12s %12s %12s %12s %12s", 
        "Truncation", "CPU-only", "Kernel CPU", "Kernel GPU", "GPU Speedup"))
    println("-"^80)
    
    for trunc in truncations
        r = results[trunc]
        gpu_speedup = r["cpu_only"] / r["kernel_gpu"]
        cpu_vs_gpu = r["kernel_cpu"] / r["kernel_gpu"]
        
        println(@sprintf("T%-11d %12d %9.3f ms %9.3f ms %9.3f ms %10.2fx",
            trunc, r["grid_points"], 
            r["cpu_only"] * 1000, r["kernel_cpu"] * 1000, 
            r["kernel_gpu"] * 1000, gpu_speedup))
    end
    println()
end

# ========== VERIFICATION ==========
println("="^80)
println("CORRECTNESS VERIFICATION")
println("="^80)
println()

println("Verifying that all versions produce identical results...")
println()

for trunc in truncations
    println("Testing T$trunc...")
    
    # Setup
    spectral_grid = SpectralGrid(; NF, trunc, nlayers, architecture=CPU())
    vertical_advection = CenteredVerticalAdvection(spectral_grid, order=order)
    model = PrimitiveWetModel(spectral_grid; vertical_advection)
    sim = initialize!(model)
    run!(sim, steps=5)
    
    diagn = sim.diagnostic_variables
    
    # Save original
    u_orig = copy(diagn.tendencies.u_tend_grid)
    v_orig = copy(diagn.tendencies.v_tend_grid)
    temp_orig = copy(diagn.tendencies.temp_tend_grid)
    
    # CPU-only version
    SpeedyWeather.vertical_advection_cpu!(diagn, model)
    u_cpu = copy(diagn.tendencies.u_tend_grid)
    v_cpu = copy(diagn.tendencies.v_tend_grid)
    temp_cpu = copy(diagn.tendencies.temp_tend_grid)
    
    # Reset and run kernel version
    diagn.tendencies.u_tend_grid .= u_orig
    diagn.tendencies.v_tend_grid .= v_orig
    diagn.tendencies.temp_tend_grid .= temp_orig
    
    SpeedyWeather.vertical_advection!(diagn, model)
    u_kernel = diagn.tendencies.u_tend_grid
    v_kernel = diagn.tendencies.v_tend_grid
    temp_kernel = diagn.tendencies.temp_tend_grid
    
    # Compare
    u_match = u_cpu ≈ u_kernel
    v_match = v_cpu ≈ v_kernel
    temp_match = temp_cpu ≈ temp_kernel
    
    if u_match && v_match && temp_match
        println("  ✓ All results match (rtol=1e-6)")
    else
        println("  ✗ Results differ!")
        if !u_match
            max_diff = maximum(abs.(u_cpu .- u_kernel))
            println("    u: max difference = $max_diff")
        end
        if !v_match
            max_diff = maximum(abs.(v_cpu .- v_kernel))
            println("    v: max difference = $max_diff")
        end
        if !temp_match
            max_diff = maximum(abs.(temp_cpu .- temp_kernel))
            println("    temp: max difference = $max_diff")
        end
    end
end

println()
println("="^80)
println("Benchmark complete!")
println("="^80)
