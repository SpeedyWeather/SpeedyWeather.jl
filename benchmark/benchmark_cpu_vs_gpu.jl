#!/usr/bin/env julia

# Benchmark script to compare CPU vs GPU performance for SpeedyWeather.jl simulations
# across different truncations

import Pkg 
Pkg.activate("benchmark")

using BenchmarkTools
using Statistics
using Printf
using CUDA
using SpeedyWeather
using Plots
using Dates

"""
    setup_simulation(trunc::Int; nlayers::Int=1, arch=CPU())

Set up a simulation with the given truncation and architecture.
Returns a simulation object ready to be run.
"""
function setup_simulation(trunc::Int; nlayers::Int=1, arch=SpeedyWeather.CPU())
    # Create a spectral grid with the given truncation and architecture
    spectral_grid = SpectralGrid(trunc=trunc, nlayers=nlayers, architecture=arch)
    
    # Create a model with the spectral grid
    model = BarotropicModel(
        spectral_grid=spectral_grid,
    )
    
    # Create a simulation with the model
    simulation = initialize!(model)
    
    # spin up 
    run!(simulation, period=Day(10))
    
    return simulation
end

"""
    run_benchmark(simulation)

Run a benchmark for the given simulation.
Returns the median time in milliseconds.
"""
function run_benchmark(simulation; ntrials=10, nsteps=100)
    
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
function benchmark_cpu_vs_gpu(truncations::Vector{Int}; nlayers::Int=1, ntrials::Int=10, nsteps::Int=100)
    cpu_times = Float64[]
    gpu_times = Float64[]
    
    println("Benchmarking CPU vs GPU performance for different truncations")
    println("===========================================================")
    println("Truncation | CPU Time (ms) | GPU Time (ms) | Speedup")
    println("----------|--------------|--------------|--------")
    
    for trunc in truncations
        # Skip GPU benchmarks if CUDA is not available
        if !CUDA.functional()
            println("CUDA is not available, skipping GPU benchmarks")
            break
        end
        
        # CPU benchmark
        cpu_sim = setup_simulation(trunc, nlayers=nlayers, arch=SpeedyWeather.CPU())
        cpu_time = run_benchmark(cpu_sim; ntrials=ntrials, nsteps=nsteps)
        push!(cpu_times, cpu_time)
        
        # GPU benchmark 
        gpu_sim = setup_simulation(trunc, nlayers=nlayers, arch=SpeedyWeather.GPU())
        gpu_time = run_benchmark(gpu_sim; ntrials=ntrials, nsteps=nsteps)
        push!(gpu_times, gpu_time)
        
        # Calculate speedup
        speedup = cpu_time / gpu_time
        
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
        ylabel="Time (ms)",
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
        ylabel="Speedup (CPU time / GPU time)",
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

# Main benchmark function
function main()
    # Define truncations to test
    truncations = [31, 63, 127, 255, 511]
    
    # Single layer benchmarks
    cpu_times, gpu_times = benchmark_cpu_vs_gpu(truncations)
    
    # Plot results
    plot_results(truncations, cpu_times, gpu_times)
end

# Run the benchmarks
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
