using SpeedyWeather
using BenchmarkTools
using Printf
using KernelAbstractions
using LinearAlgebra

"""
Benchmark script to compare the performance of horizontal_diffusion! (loop-based)
and horizontal_diffusion_kernel! (kernel-based) implementations on CPU.
"""

function create_diffusion_benchmark(; truncs=[31, 63, 127, 255], nlayers=[1, 10, 20, 40])
    suite = BenchmarkGroup()
    
    # Create benchmark group for CPU
    suite["CPU"] = BenchmarkGroup()
    
    for trunc in truncs
        suite["CPU"][string(trunc)] = BenchmarkGroup()
        
        for nlayer in nlayers
            # Create spectral grid and transform
            spectral_grid = SpectralGrid(trunc=trunc, nlayers=nlayer)
            
            # Create a model and initialize diffusion to get proper matrices
            if nlayer == 1
                # Use barotropic model for single layer
                model = BarotropicModel(spectral_grid=spectral_grid)
                diffusion = HyperDiffusion(spectral_grid)
                initialize!(diffusion, model)
            else
                # Use primitive wet model for multiple layers
                model = PrimitiveWetModel(spectral_grid=spectral_grid)
                diffusion = HyperDiffusion(spectral_grid)
                initialize!(diffusion, model)
            end
            
            # Get the properly initialized diffusion matrices
            expl = diffusion.expl
            impl = diffusion.impl
            
            # Create random spectral fields
            var = rand(Float64, spectral_grid.spectrum, nlayer)
            tendency = rand(Float64, spectral_grid.spectrum, nlayer)
            
            # Create benchmarks for CPU
            suite["CPU"][string(trunc)][string(nlayer)] = BenchmarkGroup()
            suite["CPU"][string(trunc)][string(nlayer)]["original"] = @benchmarkable SpeedyWeather.horizontal_diffusion!(
                $(copy(tendency)), $(copy(var)), $expl, $impl
            )
            suite["CPU"][string(trunc)][string(nlayer)]["kernel"] = @benchmarkable SpeedyWeather.horizontal_diffusion_kernel!(
                $(copy(tendency)), $(copy(var)), $expl, $impl
            )
        end
    end
    
    return suite
end

function run_diffusion_benchmarks()
    # Create benchmark suite
    suite = create_diffusion_benchmark()
    
    # Tune and run benchmarks
    tune!(suite)
    results = run(suite, verbose=true)
    
    # Print results in a formatted table
    println("\n\n========== HORIZONTAL DIFFUSION BENCHMARK RESULTS ==========")
    println("Format: median time in milliseconds")
    println("\n--- CPU RESULTS ---")
    
    # Print CPU results
    for trunc in sort([parse(Int, k) for k in keys(results["CPU"])])
        println("\nTruncation = $trunc")
        println("Layers | Original Implementation | Kernel Implementation | Speedup")
        println("------|------------------------|---------------------|--------")
        
        for nlayer in sort([parse(Int, k) for k in keys(results["CPU"][string(trunc)])])
            orig_time = median(results["CPU"][string(trunc)][string(nlayer)]["original"]).time / 1e6
            kernel_time = median(results["CPU"][string(trunc)][string(nlayer)]["kernel"]).time / 1e6
            speedup = orig_time / kernel_time
            
            @printf("%6d | %12.3f ms | %12.3f ms | %6.2fx\n", 
                    nlayer, orig_time, kernel_time, speedup)
        end
    end
    
    # Print summary statistics
    println("\n--- SUMMARY ---")
    println("Average speedup across all configurations:")
    
    total_speedup = 0.0
    count = 0
    
    for trunc in keys(results["CPU"])
        for nlayer in keys(results["CPU"][trunc])
            orig_time = median(results["CPU"][trunc][nlayer]["original"]).time / 1e6
            kernel_time = median(results["CPU"][trunc][nlayer]["kernel"]).time / 1e6
            speedup = orig_time / kernel_time
            total_speedup += speedup
            count += 1
        end
    end
    
    avg_speedup = total_speedup / count
    @printf("Average speedup: %.2fx\n", avg_speedup)
    
    # Find best and worst cases
    best_speedup = 0.0
    best_config = ("", "")
    worst_speedup = Inf
    worst_config = ("", "")
    
    for trunc in keys(results["CPU"])
        for nlayer in keys(results["CPU"][trunc])
            orig_time = median(results["CPU"][trunc][nlayer]["original"]).time / 1e6
            kernel_time = median(results["CPU"][trunc][nlayer]["kernel"]).time / 1e6
            speedup = orig_time / kernel_time
            
            if speedup > best_speedup
                best_speedup = speedup
                best_config = (trunc, nlayer)
            end
            
            if speedup < worst_speedup
                worst_speedup = speedup
                worst_config = (trunc, nlayer)
            end
        end
    end
    
    @printf("Best speedup: %.2fx (trunc=%s, layers=%s)\n", 
            best_speedup, best_config[1], best_config[2])
    @printf("Worst speedup: %.2fx (trunc=%s, layers=%s)\n", 
            worst_speedup, worst_config[1], worst_config[2])
    
    # Check for correctness
    println("\n--- CORRECTNESS CHECK ---")
    println("Checking if both implementations produce the same results...")
    
    # Use a moderate size for correctness check
    spectral_grid = SpectralGrid(trunc=63, nlayers=10)
    
    # Create a model and initialize diffusion to get proper matrices
    model = PrimitiveWetModel(spectral_grid=spectral_grid)
    diffusion = HyperDiffusion(spectral_grid)
    initialize!(diffusion, model)
    
    # Get the properly initialized diffusion matrices
    expl = diffusion.expl
    impl = diffusion.impl
    
    # Create random spectral fields
    var = rand(Float64, spectral_grid.spectrum, 10)
    tendency_orig = rand(Float64, spectral_grid.spectrum, 10)
    tendency_kernel = copy(tendency_orig)
    
    # Apply both implementations
    SpeedyWeather.Dynamics.horizontal_diffusion!(tendency_orig, var, expl, impl)
    SpeedyWeather.Dynamics.horizontal_diffusion_kernel!(tendency_kernel, var, expl, impl)
    
    # Check if results are identical
    max_diff = maximum(abs.(tendency_orig - tendency_kernel))
    is_identical = isapprox(tendency_orig, tendency_kernel, rtol=1e-10)
    
    if is_identical
        println("✓ Both implementations produce identical results (max difference: $max_diff)")
    else
        println("✗ Implementations produce different results (max difference: $max_diff)")
    end
end

# Run the benchmarks
println("Running horizontal diffusion benchmarks...")
run_diffusion_benchmarks()
println("All benchmarks completed.")
