"""
Benchmark script to compare particle advection implementations:
1. CPU version with for loops (particle_advection_cpu!)
2. Kernel version on CPU (particle_advection!)
3. Kernel version on GPU (particle_advection!) - if GPU available

Tests both correctness and performance.
"""

using SpeedyWeather
using BenchmarkTools
using CUDA
using Printf

# Check if GPU is available
gpu_available = CUDA.functional()

if gpu_available
    println("✓ GPU (CUDA) detected and functional")
else
    println("✗ CUDA installed but not functional")
end

println("\n" * "="^70)
println("PARTICLE ADVECTION BENCHMARK")
println("="^70)

# Test configurations
test_configs = [
    (nparticles=1000, nlayers=1, trunc=31, Grid=OctahedralGaussianGrid),
    (nparticles=10000, nlayers=1, trunc=42, Grid=OctahedralGaussianGrid),
    (nparticles=100000, nlayers=1, trunc=63, Grid=OctahedralGaussianGrid),
]

"""Helper function to setup model and run particle advection"""
function setup_and_run(nparticles, nlayers, trunc, Grid, architecture, use_cpu_version=false)
    # Create spectral grid
    spectral_grid = SpectralGrid(;
        trunc,
        nlayers,
        Grid,
        nparticles,
        architecture,
    )
    
    # Create a simple barotropic model
    model = PrimitiveWetModel(; spectral_grid, particle_advection=ParticleAdvection2D(spectral_grid, layer=1))
    
    # Initialize
    simulation = initialize!(model)
    
    # Get particles and diagnostic variables
    particles = simulation.prognostic_variables.particles
    diagn = simulation.diagnostic_variables
    clock = simulation.prognostic_variables.clock
    particle_advection = model.particle_advection
    
    # Set clock to trigger particle advection
    clock.timestep_counter = particle_advection.every_n_timesteps - 1
    
    # Make sure some particles are active
    for i in eachindex(particles)
        particles[i] = SpeedyWeather.activate(particles[i])
    end

    # Set initial particle positions
    simulation.progn.particles .= rand(Particle{Float32}, nparticles)
    
    # Choose which version to run
    if use_cpu_version
        # Use CPU for-loop version
        return () -> SpeedyWeather.particle_advection_cpu!(
            particles, diagn, clock, particle_advection
        ), particles, diagn
    else
        # Use kernel version
        return () -> SpeedyWeather.particle_advection!(
            particles, diagn, clock, particle_advection
        ), particles, diagn
    end
end

"""Compare correctness of two implementations"""
function compare_correctness(nparticles, nlayers, trunc, Grid)
    println("\n" * "-"^70)
    println("Correctness Test: $nparticles particles, T$trunc")
    println("-"^70)
    
    # Setup CPU for-loop version
    run_cpu!, particles_cpu, diagn_cpu = setup_and_run(
        nparticles, nlayers, trunc, Grid, CPU(), true
    )
    
    # Setup CPU kernel version
    run_kernel_cpu!, particles_kernel, diagn_kernel = setup_and_run(
        nparticles, nlayers, trunc, Grid, CPU(), false
    )
    
    # Copy initial state to both
    for i in eachindex(particles_cpu, particles_kernel)
        particles_kernel[i] = particles_cpu[i]
    end
    
    # Run both versions
    run_cpu!()
    run_kernel_cpu!()
    
    # Compare results
    max_diff_lon = 0.0
    max_diff_lat = 0.0
    n_compared = 0
    
    for i in eachindex(particles_cpu, particles_kernel)
        if SpeedyWeather.isactive(particles_cpu[i]) && SpeedyWeather.isactive(particles_kernel[i])
            diff_lon = abs(particles_cpu[i].lon - particles_kernel[i].lon)
            diff_lat = abs(particles_cpu[i].lat - particles_kernel[i].lat)
            max_diff_lon = max(max_diff_lon, diff_lon)
            max_diff_lat = max(max_diff_lat, diff_lat)
            n_compared += 1
        end
    end
    
    println("Compared $n_compared active particles")
    println("Max difference in longitude: $(@sprintf("%.2e", max_diff_lon))°")
    println("Max difference in latitude:  $(@sprintf("%.2e", max_diff_lat))°")
    
    # Check if results match within tolerance
    tolerance = 1e-10
    if max_diff_lon < tolerance && max_diff_lat < tolerance
        println("✓ Results match within tolerance ($(@sprintf("%.2e", tolerance))°)")
        return true
    else
        println("✗ Results differ beyond tolerance!")
        return false
    end
end

"""Benchmark a single configuration"""
function benchmark_config(nparticles, nlayers, trunc, Grid)
    println("\n" * "="^70)
    println("Benchmarking: $nparticles particles, T$trunc, $(Grid)")
    println("="^70)
    
    # Test correctness first
    correct = compare_correctness(nparticles, nlayers, trunc, Grid)
    if !correct
        println("⚠ Warning: Correctness test failed, but continuing with benchmarks...")
    end
    
    println("\n" * "-"^70)
    println("Performance Benchmarks")
    println("-"^70)
    
    # Benchmark CPU for-loop version
    println("\n1. CPU For-Loop Version:")
    run_cpu!, _, _ = setup_and_run(nparticles, nlayers, trunc, Grid, CPU(), true)
    bench_cpu = @benchmark $run_cpu!() samples=100 seconds=10
    println(bench_cpu)
    time_cpu = median(bench_cpu).time / 1e6  # Convert to ms
    
    # Benchmark CPU kernel version
    println("\n2. CPU Kernel Version:")
    run_kernel_cpu!, _, _ = setup_and_run(nparticles, nlayers, trunc, Grid, CPU(), false)
    bench_kernel_cpu = @benchmark $run_kernel_cpu!() samples=100 seconds=10
    println(bench_kernel_cpu)
    time_kernel_cpu = median(bench_kernel_cpu).time / 1e6  # Convert to ms
    
    # Benchmark GPU kernel version if available
    time_kernel_gpu = nothing
    if gpu_available
        println("\n3. GPU Kernel Version:")
        try
            run_kernel_gpu!, _, _ = setup_and_run(nparticles, nlayers, trunc, Grid, GPU(), false)
            # Warmup
            run_kernel_gpu!()
            CUDA.synchronize()
            bench_kernel_gpu = @benchmark begin
                $run_kernel_gpu!()
                CUDA.synchronize()
            end samples=100 seconds=10
            println(bench_kernel_gpu)
            time_kernel_gpu = median(bench_kernel_gpu).time / 1e6  # Convert to ms
        catch e
            println("✗ GPU benchmark failed: $e")
        end
    end
    
    # Summary
    println("\n" * "-"^70)
    println("Summary for $nparticles particles:")
    println("-"^70)
    println(@sprintf("CPU For-Loop:    %8.3f ms", time_cpu))
    println(@sprintf("CPU Kernel:      %8.3f ms  (%.2fx vs for-loop)", 
        time_kernel_cpu, time_cpu / time_kernel_cpu))
    
    if time_kernel_gpu !== nothing
        println(@sprintf("GPU Kernel:      %8.3f ms  (%.2fx vs for-loop, %.2fx vs CPU kernel)", 
            time_kernel_gpu, time_cpu / time_kernel_gpu, time_kernel_cpu / time_kernel_gpu))
    end
    
    return (
        cpu_forloop = time_cpu,
        cpu_kernel = time_kernel_cpu,
        gpu_kernel = time_kernel_gpu,
        speedup_kernel_cpu = time_cpu / time_kernel_cpu,
        speedup_gpu = time_kernel_gpu !== nothing ? time_cpu / time_kernel_gpu : nothing
    )
end

# Run benchmarks for all configurations
println("\nStarting benchmarks...")
results = []

for config in test_configs
    result = benchmark_config(config.nparticles, config.nlayers, config.trunc, config.Grid)
    push!(results, (nparticles=config.nparticles, result...))
end

# Final summary
println("\n" * "="^70)
println("FINAL SUMMARY")
println("="^70)
println("\n" * @sprintf("%-12s | %10s | %10s | %10s | %8s", 
    "Particles", "CPU Loop", "CPU Kernel", "GPU Kernel", "Speedup"))
println("-"^70)

for r in results
    gpu_str = r.gpu_kernel !== nothing ? @sprintf("%.3f ms", r.gpu_kernel) : "N/A"
    speedup_str = r.speedup_gpu !== nothing ? @sprintf("%.2fx", r.speedup_gpu) : "N/A"
    
    println(@sprintf("%-12d | %8.3f ms | %8.3f ms | %10s | %8s",
        r.nparticles, r.cpu_forloop, r.cpu_kernel, gpu_str, speedup_str))
end

println("\n" * "="^70)
println("Benchmark complete!")
println("="^70)
