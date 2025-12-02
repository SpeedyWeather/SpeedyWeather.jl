"""
Benchmark and Correctness Test for Kernelized Tendency Functions

Simplified version focusing on PrimitiveWetModel only.
Tests performance of kernelized vs original implementations on CPU.
"""
import Pkg 
Pkg.activate(".")

using SpeedyWeather
using BenchmarkTools
using Printf

# Create test model
function create_test_model(; trunc=31, nlayers=8)
    spectral_grid = SpectralGrid(; trunc, nlayers, Grid=FullGaussianGrid, NF=Float64)
    model = PrimitiveWetModel(; spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps=4) # small spinup 
    return simulation
end

# Benchmark a function
function benchmark_function(func!, args...; samples=50, evals=5)
    # Warmup
    for _ in 1:2
        func!(args...)
    end
    
    # Benchmark
    b = @benchmark $func!($(args)...) samples=samples evals=evals
    return b
end

# Check correctness by comparing outputs
function check_correctness(name::String, func_original!, func_kernel!, args_orig, args_kern; rtol=1e-10, atol=1e-12)
    println("\n  Checking correctness...")
    
    # Run both versions
    func_original!(args_orig...)
    func_kernel!(args_kern...)
    
    # Compare the modified arrays (assuming they're in diagn)
    # This is a simplified check - we compare key output arrays
    diagn_orig = args_orig[1]
    diagn_kern = args_kern[1]
    
    # Check if results are approximately equal
    all_close = true
    max_diff = 0.0
    
    # Compare a few key fields that are likely modified
    for field in [:u_grid, :v_grid, :div_grid, :temp_grid]
        if hasfield(typeof(diagn_orig.grid), field)
            arr_orig = getfield(diagn_orig.grid, field)
            arr_kern = getfield(diagn_kern.grid, field)
            
            diff = maximum(abs.(arr_orig .- arr_kern))
            max_diff = max(max_diff, diff)
            
            if !isapprox(arr_orig, arr_kern; rtol=rtol, atol=atol)
                all_close = false
            end
        end
    end
    
    if all_close
        println("  ✓ Correctness check PASSED (max diff: $(max_diff))")
    else
        println("  ✗ Correctness check FAILED (max diff: $(max_diff))")
    end
    
    return all_close
end

# Compare performance
function compare_performance(name::String, func_original!, func_kernel!, args...)
    println("\n" * "="^70)
    println("Benchmarking: $name")
    println("="^70)
    
    # First check correctness with fresh copies of the state
    # We need to create separate copies for each function
    # For now, we'll skip detailed correctness checking and just run both
    
    b_original = benchmark_function(func_original!, args...)
    b_kernel = benchmark_function(func_kernel!, args...)
    
    # Display results
    t_orig = median(b_original).time
    t_kern = median(b_kernel).time
    speedup = t_orig / t_kern
    
    println("  Original:   $(BenchmarkTools.prettytime(t_orig))")
    println("  Kernelized: $(BenchmarkTools.prettytime(t_kern))")
    println("  Speedup:    $(round(speedup, digits=2))x")
    
    if speedup > 1.10
        println("  ✓ Kernelized version is faster!")
    elseif speedup > 0.90
        println("  ≈ Performance is similar")
    else
        println("  ⚠ Original version is faster")
    end
    
    return (original=b_original, kernel=b_kernel, speedup=speedup)
end

# Test correctness for a specific function
function test_function_correctness(name::String, func_original!, func_kernel!, trunc, nlayers)
    println("\n" * "-"^70)
    println("Testing correctness: $name")
    println("-"^70)
    
    # Create two separate simulations
    sim_orig = create_test_model(; trunc, nlayers)
    sim_kern = create_test_model(; trunc, nlayers)
    
    progn_orig, diagn_orig, model_orig = SpeedyWeather.unpack(sim_orig)
    progn_kern, diagn_kern, model_kern = SpeedyWeather.unpack(sim_kern)
    lf = 1
    
    # Build argument lists for each function
    args_orig, args_kern = if name == "vertical_integration!"
        ([diagn_orig, progn_orig, lf, model_orig.geometry],
         [diagn_kern, progn_kern, lf, model_kern.geometry])
    elseif name == "vertical_velocity!"
        ([diagn_orig, model_orig.geometry],
         [diagn_kern, model_kern.geometry])
    elseif name == "vordiv_tendencies!"
        ([diagn_orig, model_orig.coriolis, model_orig.atmosphere, model_orig.geometry, model_orig.spectral_transform],
         [diagn_kern, model_kern.coriolis, model_kern.atmosphere, model_kern.geometry, model_kern.spectral_transform])
    elseif name == "temperature_tendency!"
        ([diagn_orig, model_orig.adiabatic_conversion, model_orig.atmosphere, model_orig.implicit, model_orig.geometry, model_orig.spectral_transform],
         [diagn_kern, model_kern.adiabatic_conversion, model_kern.atmosphere, model_kern.implicit, model_kern.geometry, model_kern.spectral_transform])
    else
        return false  # Skip for now
    end
    
    # Run both versions
    func_original!(args_orig...)
    func_kernel!(args_kern...)
    
    # Compare results - check all arrays in diagn.grid, diagn.dynamics, and diagn.tendencies
    max_diff = 0.0
    all_close = true
    failed_fields = String[]
    
    # Helper function to compare fields recursively
    function compare_fields(obj_orig, obj_kern, prefix="")
        for field in fieldnames(typeof(obj_orig))
            arr_orig = getfield(obj_orig, field)
            arr_kern = getfield(obj_kern, field)
            field_name = prefix * string(field)
            
            if arr_orig isa AbstractArray && arr_kern isa AbstractArray
                diff = maximum(abs.(arr_orig .- arr_kern))
                max_diff = max(max_diff, diff)
                
                if !isapprox(arr_orig, arr_kern; rtol=1e-10, atol=1e-12)
                    all_close = false
                    push!(failed_fields, "$field_name (max diff: $diff)")
                end
            end
        end
    end
    
    # Compare grid fields
    compare_fields(diagn_orig.grid, diagn_kern.grid, "grid.")
    
    # Compare dynamics fields
    compare_fields(diagn_orig.dynamics, diagn_kern.dynamics, "dynamics.")
    
    # Compare tendencies fields
    compare_fields(diagn_orig.tendencies, diagn_kern.tendencies, "tendencies.")
    
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

# Main benchmark suite
function run_benchmarks(; trunc=31, nlayers=8, check_correctness=true)
    println("\n" * "="^70)
    println("KERNELIZED TENDENCIES BENCHMARK")
    println("="^70)
    println("Model: PrimitiveWetModel")
    println("Truncation: T$trunc")
    println("Layers: $nlayers")
    println("Architecture: CPU")
    
    # Run correctness tests first if requested
    if check_correctness
        println("\n" * "="^70)
        println("CORRECTNESS TESTS")
        println("="^70)
        
        correctness_results = Dict{String, Bool}()
        correctness_results["vertical_integration!"] = test_function_correctness(
            "vertical_integration!", 
            SpeedyWeather.vertical_integration!,
            SpeedyWeather.vertical_integration_kernel!,
            trunc, nlayers
        )
        correctness_results["vertical_velocity!"] = test_function_correctness(
            "vertical_velocity!",
            SpeedyWeather.vertical_velocity!,
            SpeedyWeather.vertical_velocity_kernel!,
            trunc, nlayers
        )
        correctness_results["vordiv_tendencies!"] = test_function_correctness(
            "vordiv_tendencies!",
            SpeedyWeather.vordiv_tendencies!,
            SpeedyWeather.vordiv_tendencies_kernel!,
            trunc, nlayers
        )
        correctness_results["temperature_tendency!"] = test_function_correctness(
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
    println("PERFORMANCE BENCHMARKS")
    println("="^70)
    
    simulation = create_test_model(; trunc, nlayers)
    progn, diagn, model = SpeedyWeather.unpack(simulation)
    lf = 1
    
    results = Dict{String, Any}()
    
    # Test 1: vertical_integration!
    results["vertical_integration"] = compare_performance(
        "vertical_integration!",
        SpeedyWeather.vertical_integration!,
        SpeedyWeather.vertical_integration_kernel!,
        diagn, progn, lf, model.geometry
    )
    
    # Test 2: vertical_velocity!
    results["vertical_velocity"] = compare_performance(
        "vertical_velocity!",
        SpeedyWeather.vertical_velocity!,
        SpeedyWeather.vertical_velocity_kernel!,
        diagn, model.geometry
    )
    
    # Test 3: vordiv_tendencies!
    results["vordiv_tendencies"] = compare_performance(
        "vordiv_tendencies!",
        SpeedyWeather.vordiv_tendencies!,
        SpeedyWeather.vordiv_tendencies_kernel!,
        diagn, model.coriolis, model.atmosphere, model.geometry, model.spectral_transform
    )
    
    # Test 4: temperature_tendency!
    results["temperature_tendency"] = compare_performance(
        "temperature_tendency!",
        SpeedyWeather.temperature_tendency!,
        SpeedyWeather.temperature_tendency_kernel!,
        diagn, model.adiabatic_conversion, model.atmosphere, model.implicit,
        model.geometry, model.spectral_transform
    )
    
    # Test 5: pressure_gradient_flux!
    results["pressure_gradient_flux"] = compare_performance(
        "pressure_gradient_flux!",
        SpeedyWeather.pressure_gradient_flux!,
        SpeedyWeather.pressure_gradient_flux_kernel!,
        diagn, progn, lf, model.spectral_transform
    )
    
    # Test 6: temperature_anomaly!
    results["temperature_anomaly"] = compare_performance(
        "temperature_anomaly!",
        SpeedyWeather.temperature_anomaly!,
        SpeedyWeather.temperature_anomaly_kernel!,
        diagn, model.implicit
    )
    
    # Test 7: horizontal_advection!
    A_tend = diagn.tendencies.temp_tend
    A_tend_grid = diagn.tendencies.temp_tend_grid
    A_grid = diagn.grid.temp_grid
    
    results["horizontal_advection"] = compare_performance(
        "horizontal_advection!",
        SpeedyWeather.horizontal_advection!,
        SpeedyWeather.horizontal_advection_kernel!,
        A_tend, A_tend_grid, A_grid, diagn, model.geometry, model.spectral_transform
    )
    
    # Test 8: flux_divergence!
    results["flux_divergence"] = compare_performance(
        "flux_divergence!",
        SpeedyWeather.flux_divergence!,
        SpeedyWeather.flux_divergence_kernel!,
        A_tend, A_grid, diagn, model.geometry, model.spectral_transform
    )
    
    # Test 9: linear_pressure_gradient!
    results["linear_pressure_gradient"] = compare_performance(
        "linear_pressure_gradient!",
        SpeedyWeather.linear_pressure_gradient!,
        SpeedyWeather.linear_pressure_gradient_kernel!,
        diagn, progn, lf, model.atmosphere, model.implicit
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
    
    avg_speedup = mean([r.speedup for r in values(results)])
    println("\nAverage speedup: $(round(avg_speedup, digits=2))x")
    
    return results
end

# Quick test with small grid
function quick_test()
    run_benchmarks(; trunc=15, nlayers=4)
end

# Full benchmark with multiple sizes
function full_benchmark()
    configs = [
        (trunc=31, nlayers=8, name="Small (T31, 8 layers)"),
        (trunc=63, nlayers=16, name="Medium (T63, 16 layers)"),
    ]
    
    all_results = Dict{String, Any}()
    
    for config in configs
        println("\n\n" * "█"^70)
        println("Configuration: $(config.name)")
        println("█"^70)
        
        results = run_benchmarks(; trunc=config.trunc, nlayers=config.nlayers)
        all_results[config.name] = results
    end
    
    # Summary across configurations
    println("\n\n" * "="^70)
    println("CROSS-CONFIGURATION SUMMARY")
    println("="^70)
    
    for (config_name, results) in all_results
        avg_speedup = mean([r.speedup for r in values(results)])
        @printf("%-30s: %.2fx average\n", config_name, avg_speedup)
    end
    
    return all_results
end

# Export main functions
export run_benchmarks, quick_test, full_benchmark

# Run quick test if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running quick benchmark test...")
    quick_test()
end
