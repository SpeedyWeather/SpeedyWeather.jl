import Pkg 
Pkg.activate("test")

using SpeedyWeather
using BenchmarkTools
using Printf
using CUDA  # For GPU support

# Check if GPU is available
const HAS_GPU = CUDA.functional()

if !HAS_GPU
    @warn "No CUDA-capable GPU detected. This script requires a GPU to run."
    exit(1)
end

# Create test model with specified architecture
function create_test_model(; trunc=42, nlayers=8, architecture=CPU())
    spectral_grid = SpectralGrid(; trunc, nlayers, NF=Float32, architecture)
    model = PrimitiveWetModel(; spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps=4) # small spinup 
    return simulation
end

# Benchmark a function
function benchmark_function(func!, args...; samples=100, evals=5)
    # Warmup
    for _ in 1:2
        func!(args...)
    end
    
    # Synchronize GPU if needed
    if any(arg -> arg isa CuArray || (hasproperty(arg, :data) && arg.data isa CuArray), args)
        CUDA.synchronize()
    end
    
    # Benchmark
    b = @benchmark begin
        $func!($(args)...)
        if any(arg -> arg isa CuArray || (hasproperty(arg, :data) && arg.data isa CuArray), $args)
            CUDA.synchronize()
        end
    end samples=samples evals=evals
    
    return b
end

# Test correctness by comparing CPU and GPU results
function test_function_correctness_gpu(name::String, func_cpu!, func_gpu!, trunc, nlayers)
    println("\n" * "-"^70)
    println("Testing correctness (CPU vs GPU): $name")
    println("-"^70)
    
    # Create two separate simulations - one on CPU, one on GPU
    sim_cpu = create_test_model(; trunc, nlayers, architecture=CPU())
    sim_gpu = create_test_model(; trunc, nlayers, architecture=GPU())
    
    progn_cpu, diagn_cpu, model_cpu = SpeedyWeather.unpack(sim_cpu)
    progn_gpu, diagn_gpu, model_gpu = SpeedyWeather.unpack(sim_gpu)
    lf = 1
    
    # Build argument lists for each function
    args_cpu, args_gpu = if name == "vertical_integration!"
        ([diagn_cpu, progn_cpu, lf, model_cpu.geometry],
         [diagn_gpu, progn_gpu, lf, model_gpu.geometry])
    elseif name == "vertical_velocity!"
        ([diagn_cpu, model_cpu.geometry],
         [diagn_gpu, model_gpu.geometry])
    elseif name == "vordiv_tendencies!"
        ([diagn_cpu, model_cpu.coriolis, model_cpu.atmosphere, model_cpu.geometry, model_cpu.spectral_transform],
         [diagn_gpu, model_gpu.coriolis, model_gpu.atmosphere, model_gpu.geometry, model_gpu.spectral_transform])
    elseif name == "temperature_tendency!"
        ([diagn_cpu, model_cpu.adiabatic_conversion, model_cpu.atmosphere, model_cpu.implicit, model_cpu.geometry, model_cpu.spectral_transform],
         [diagn_gpu, model_gpu.adiabatic_conversion, model_gpu.atmosphere, model_gpu.implicit, model_gpu.geometry, model_gpu.spectral_transform])
    else
        return false  # Skip for now
    end
    
    # Run both versions
    func_cpu!(args_cpu...)
    func_gpu!(args_gpu...)
    CUDA.synchronize()
    
    # Compare results - check all arrays in diagn.grid, diagn.dynamics, and diagn.tendencies
    max_diff = 0.0
    all_close = true
    failed_fields = String[]
    
    # Helper function to compare fields recursively
    function compare_fields(obj_cpu, obj_gpu, prefix="")
        for field in fieldnames(typeof(obj_cpu))
            arr_cpu = getfield(obj_cpu, field)
            arr_gpu = getfield(obj_gpu, field)
            field_name = prefix * string(field)
            
            if arr_cpu isa AbstractArray && arr_gpu isa AbstractArray
                # Transfer GPU array to CPU for comparison
                arr_gpu_cpu = on_architecture(CPU(), arr_gpu)
                
                diff = maximum(abs.(arr_cpu .- arr_gpu_cpu))
                max_diff = max(max_diff, diff)
                
                # Use slightly relaxed tolerance for GPU (floating point differences)
                if !isapprox(arr_cpu, arr_gpu_cpu; rtol=1e-8, atol=1e-10)
                    all_close = false
                    push!(failed_fields, "$field_name (max diff: $diff)")
                end
            end
        end
    end
    
    # Compare grid fields
    compare_fields(diagn_cpu.grid, diagn_gpu.grid, "grid.")
    
    # Compare dynamics fields
    compare_fields(diagn_cpu.dynamics, diagn_gpu.dynamics, "dynamics.")
    
    # Compare tendencies fields
    compare_fields(diagn_cpu.tendencies, diagn_gpu.tendencies, "tendencies.")
    
    if all_close
        println("  ✓ Correctness PASSED (max diff: $(max_diff))")
    else
        println("  ✗ Correctness FAILED (max diff: $(max_diff))")
        println("  Failed fields:")
        for field in failed_fields
            println("    - $field")
        end
    end
    
    return all_close
end

# Compare performance between CPU (original) and GPU (kernelized)
function compare_performance_gpu(name::String, func_cpu!, func_gpu!, args_cpu, args_gpu)
    println("\n" * "="^70)
    println("Benchmarking: $name")
    println("="^70)
    
    b_cpu = benchmark_function(func_cpu!, args_cpu...)
    b_gpu = benchmark_function(func_gpu!, args_gpu...)
    
    # Display results
    t_cpu = median(b_cpu).time
    t_gpu = median(b_gpu).time
    speedup = t_cpu / t_gpu
    
    println("  CPU (original):   $(BenchmarkTools.prettytime(t_cpu))")
    println("  GPU (kernelized): $(BenchmarkTools.prettytime(t_gpu))")
    println("  Speedup:          $(round(speedup, digits=2))x")
    
    if speedup > 1.10
        println("  ✓ GPU version is faster!")
    elseif speedup > 0.90
        println("  ≈ Performance is similar")
    else
        println("  ⚠ CPU version is faster")
    end
    
    return (cpu=b_cpu, gpu=b_gpu, speedup=speedup)
end

# Main benchmark suite
function run_benchmarks_gpu(; trunc=31, nlayers=8, check_correctness=false)
    println("\n" * "="^70)
    println("KERNELIZED TENDENCIES BENCHMARK - GPU vs CPU")
    println("="^70)
    println("Model: PrimitiveWetModel")
    println("Truncation: T$trunc")
    println("Layers: $nlayers")
    println("CPU Architecture: $(CPU())")
    println("GPU Architecture: $(GPU())")
    println("GPU Device: $(CUDA.device())")
    
    # Run correctness tests first if requested
    if check_correctness
        println("\n" * "="^70)
        println("CORRECTNESS TESTS (CPU vs GPU)")
        println("="^70)
        
        correctness_results = Dict{String, Bool}()
        correctness_results["vertical_integration!"] = test_function_correctness_gpu(
            "vertical_integration!", 
            SpeedyWeather.vertical_integration!,
            SpeedyWeather.vertical_integration_kernel!,
            trunc, nlayers
        )
        correctness_results["vertical_velocity!"] = test_function_correctness_gpu(
            "vertical_velocity!",
            SpeedyWeather.vertical_velocity!,
            SpeedyWeather.vertical_velocity_kernel!,
            trunc, nlayers
        )
        correctness_results["vordiv_tendencies!"] = test_function_correctness_gpu(
            "vordiv_tendencies!",
            SpeedyWeather.vordiv_tendencies!,
            SpeedyWeather.vordiv_tendencies_kernel!,
            trunc, nlayers
        )
        correctness_results["temperature_tendency!"] = test_function_correctness_gpu(
            "temperature_tendency!",
            SpeedyWeather.temperature_tendency!,
            SpeedyWeather.temperature_tendency_kernel!,
            trunc, nlayers
        )
        
        # Summary
        println("\n" * "-"^70)
        all_passed = all(values(correctness_results))
        if all_passed
            println("✓ All correctness tests PASSED")
        else
            println("✗ Some correctness tests FAILED:")
            for (name, passed) in correctness_results
                if !passed
                    println("  - $name")
                end
            end
        end
    end
    
    # Now run performance benchmarks
    println("\n" * "="^70)
    println("PERFORMANCE BENCHMARKS (CPU vs GPU)")
    println("="^70)
    
    sim_cpu = create_test_model(; trunc, nlayers, architecture=CPU())
    sim_gpu = create_test_model(; trunc, nlayers, architecture=GPU())
    
    progn_cpu, diagn_cpu, model_cpu = SpeedyWeather.unpack(sim_cpu)
    progn_gpu, diagn_gpu, model_gpu = SpeedyWeather.unpack(sim_gpu)
    lf = 1
    
    results = Dict{String, Any}()
    
    # Test 1: vertical_integration!
    results["vertical_integration"] = compare_performance_gpu(
        "vertical_integration!",
        SpeedyWeather.vertical_integration!,
        SpeedyWeather.vertical_integration!,
        [diagn_cpu, progn_cpu, lf, model_cpu.geometry],
        [diagn_gpu, progn_gpu, lf, model_gpu.geometry]
    )
    
    # Test 2: vertical_velocity!
    results["vertical_velocity"] = compare_performance_gpu(
        "vertical_velocity!",
        SpeedyWeather.vertical_velocity!,
        SpeedyWeather.vertical_velocity!,
        [diagn_cpu, model_cpu.geometry],
        [diagn_gpu, model_gpu.geometry]
    )
    
    # Test 3: vordiv_tendencies!
    results["vordiv_tendencies"] = compare_performance_gpu(
        "vordiv_tendencies!",
        SpeedyWeather.vordiv_tendencies!,
        SpeedyWeather.vordiv_tendencies!,
        [diagn_cpu, model_cpu],
        [diagn_gpu, model_gpu]
    )
    
    # Test 4: temperature_tendency!
    results["temperature_tendency"] = compare_performance_gpu(
        "temperature_tendency!",
        SpeedyWeather.temperature_tendency!,
        SpeedyWeather.temperature_tendency!,
        [diagn_cpu, model_cpu.adiabatic_conversion, model_cpu.atmosphere, model_cpu.implicit, model_cpu.geometry, model_cpu.spectral_transform],
        [diagn_gpu, model_gpu.adiabatic_conversion, model_gpu.atmosphere, model_gpu.implicit, model_gpu.geometry, model_gpu.spectral_transform]
    )
    
    # Test 5: pressure_gradient_flux!
    results["pressure_gradient_flux"] = compare_performance_gpu(
        "pressure_gradient_flux!",
        SpeedyWeather.pressure_gradient_flux!,
        SpeedyWeather.pressure_gradient_flux!,
        [diagn_cpu, progn_cpu, lf, model_cpu.spectral_transform],
        [diagn_gpu, progn_gpu, lf, model_gpu.spectral_transform]
    )
    
    # Test 7: horizontal_advection!
    A_tend_cpu = diagn_cpu.tendencies.temp_tend
    A_tend_grid_cpu = diagn_cpu.tendencies.temp_tend_grid
    A_grid_cpu = diagn_cpu.grid.temp_grid
    
    A_tend_gpu = diagn_gpu.tendencies.temp_tend
    A_tend_grid_gpu = diagn_gpu.tendencies.temp_tend_grid
    A_grid_gpu = diagn_gpu.grid.temp_grid
    
    results["horizontal_advection"] = compare_performance_gpu(
        "horizontal_advection!",
        SpeedyWeather.horizontal_advection!,
        SpeedyWeather.horizontal_advection!,
        [A_tend_cpu, A_tend_grid_cpu, A_grid_cpu, diagn_cpu, model_cpu.geometry, model_cpu.spectral_transform],
        [A_tend_gpu, A_tend_grid_gpu, A_grid_gpu, diagn_gpu, model_gpu.geometry, model_gpu.spectral_transform]
    )
    
    # Test 8: flux_divergence!
    results["flux_divergence"] = compare_performance_gpu(
        "flux_divergence!",
        SpeedyWeather.flux_divergence!,
        SpeedyWeather.flux_divergence!,
        [A_tend_cpu, A_grid_cpu, diagn_cpu, model_cpu.geometry, model_cpu.spectral_transform],
        [A_tend_gpu, A_grid_gpu, diagn_gpu, model_gpu.geometry, model_gpu.spectral_transform]
    )
    
    # Test 9: linear_pressure_gradient!
    results["linear_pressure_gradient"] = compare_performance_gpu(
        "linear_pressure_gradient!",
        SpeedyWeather.linear_pressure_gradient!,
        SpeedyWeather.linear_pressure_gradient!,
        [diagn_cpu, progn_cpu, lf, model_cpu.atmosphere, model_cpu.implicit],
        [diagn_gpu, progn_gpu, lf, model_gpu.atmosphere, model_gpu.implicit]
    )
    
    # Summary
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    
    for (name, result) in collect(results)
        speedup = result.speedup
        status = speedup > 1.10 ? "✓" : speedup > 0.90 ? "≈" : "⚠"
        @printf("%-30s: %6.2fx  %s\n", name, speedup, status)
    end
    
    avg_speedup = sum([r.speedup for r in values(results)]) / length(results)
    println("\nAverage speedup: $(round(avg_speedup, digits=2))x")
    
    return results
end

# Quick test with small grid
function quick_test_gpu()
    run_benchmarks_gpu(; trunc=15, nlayers=4)
end

# Full benchmark with multiple sizes
function full_benchmark_gpu()
    configs = [
        (trunc=42, nlayers=8, name="Small (T42, 8 layers)"),
        (trunc=63, nlayers=16, name="Medium (T63, 16 layers)"),
        (trunc=125, nlayers=16, name="Large (T125, 16 layers)"),
    ]
    
    all_results = Dict{String, Any}()
    
    for config in configs
        println("\n\n" * "█"^70)
        println("Configuration: $(config.name)")
        println("█"^70)
        
        results = run_benchmarks_gpu(; trunc=config.trunc, nlayers=config.nlayers)
        all_results[config.name] = results
    end
    
    # Summary across configurations
    println("\n\n" * "="^70)
    println("CROSS-CONFIGURATION SUMMARY")
    println("="^70)
    
    for (config_name, results) in all_results
        avg_speedup = sum([r.speedup for r in values(results)]) / length(results)
        @printf("%-30s: %.2fx average\n", config_name, avg_speedup)
    end
    
    return all_results
end

# Export main functions
export run_benchmarks_gpu, quick_test_gpu, full_benchmark_gpu

# Run quick test if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running quick GPU benchmark test...")
    quick_test_gpu()
end
