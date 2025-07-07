

using BenchmarkTools
using Statistics
using Printf
using CUDA
using SpeedyWeather
using Plots
using Dates

"""
    run_benchmark(simulation)

Run a benchmark for the given simulation.
Returns the median time in milliseconds.
"""
function run_benchmark(simulation; ntrials=1, nsteps=1500)
    
    sypd = zeros(Float64, ntrials)
    for i in 1:ntrials
        run!(simulation, steps=nsteps)
        time_elapsed = simulation.model.feedback.progress_meter.tlast - simulation.model.feedback.progress_meter.tinit
        sypd[i] = simulation.model.time_stepping.Δt_sec*nsteps / (time_elapsed * 365.25)
    end 

    # Return the median time of SYPD
    return median(sypd)
end

"""
    benchmark_cpu_vs_gpu(truncations::Vector{Int}; nlayers::Int=1)

Benchmark CPU vs GPU performance for the given truncations.
Returns a tuple of (cpu_times, gpu_times) in milliseconds.
"""
function benchmark_cpu_vs_gpu(truncations::Vector{Int}; nlayers::Int=1, ntrials::Int=1, nsteps::Union{Int, Vector{Int}}=1000)
    cpu_times = Float64[]
    gpu_times = Float64[]
    
    if typeof(nsteps) <: Int
        nsteps = fill(nsteps, length(truncations))
    end
    
    println("Benchmarking CPU vs GPU performance for different truncations")
    println("===========================================================")
    println("Truncation | CPU Time (SYPD) | GPU Time (SYPD) | Speedup")
    println("----------|--------------|--------------|--------")
    
    for (i, trunc) in enumerate(truncations)
        # Skip GPU benchmarks if CUDA is not available
        if !CUDA.functional()
            println("CUDA is not available, skipping GPU benchmarks")
            break
        end
        
        # CPU benchmark
        cpu_sim = setup_simulation(trunc, nlayers=nlayers, arch=SpeedyWeather.CPU())
        cpu_time = run_benchmark(cpu_sim; ntrials=ntrials, nsteps=nsteps[i])
        push!(cpu_times, cpu_time)
        
        # GPU benchmark 
        gpu_sim = setup_simulation(trunc, nlayers=nlayers, arch=SpeedyWeather.GPU())
        gpu_time = run_benchmark(gpu_sim; ntrials=ntrials, nsteps=nsteps[i])
        push!(gpu_times, gpu_time)
        
        # Calculate speedup
        speedup = gpu_time / cpu_time
        
        # Print results
        @printf("%10d | %12.2f | %12.2f | %6.2fx\n", trunc, cpu_time, gpu_time, speedup)
    end
    
    return (cpu_times, gpu_times)
end

"""    
    plot_results(truncations, cpu_times, gpu_times, filename="cpu_vs_gpu_benchmark.pdf")

Plot the benchmark results and save to a PDF file.
"""
function plot_results(truncations, cpu_times, gpu_times; filename="cpu_vs_gpu_benchmark.pdf")
    # Create a plot
    p1 = plot(
        truncations, 
        [cpu_times gpu_times], 
        label=["CPU" "GPU"],
        xlabel="Truncation",
        ylabel="Time (SYPD)",
        title="CPU vs GPU Performance",
        marker=:circle,
        markersize=4,
        linewidth=2,
        legend=:topleft,
        yscale=:log10,
        xscale=:log2,
        grid=true
    )
    
    # Calculate speedup
    speedups = cpu_times ./ gpu_times
    
    # Create a speedup plot
    p2 = plot(
        truncations, 
        speedups, 
        label="GPU Speedup",
        xlabel="Truncation",
        ylabel="Speedup (GPU/CPU SYPD)",
        title="GPU Speedup Factor",
        marker=:circle,
        markersize=4,
        linewidth=2,
        legend=:topleft,
        xscale=:log2,
        grid=true,
        color=:green
    )
    
    # Combine plots
    p = plot(p1, p2, layout=(2,1), size=(800, 800))
    
    # Add timestamp
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    annotate!(p[2], [(truncations[1], minimum(speedups), ("Generated: $timestamp", 8, :left, :bottom))])
    
    # Save to PDF
    savefig(p, filename)
    println("Results saved to $filename")
    
    return p
end

