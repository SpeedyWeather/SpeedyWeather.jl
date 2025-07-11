#!/usr/bin/env julia

# Benchmark script to compare CPU vs GPU performance for SpeedyWeather.jl simulations
# across different truncations

import Pkg 
Pkg.activate("benchmark")

include("benchmark_models.jl")

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
    run!(simulation, period=Day(3))
    
    return simulation
end

# Main benchmark function
function main()
    # Define truncations to test
    truncations = [31, 63, 127, 255, 511, 1023]
    nsteps = [1500, 1000, 500, 250, 100, 50]
    # Single layer benchmarks
    cpu_times, gpu_times = benchmark_cpu_vs_gpu(truncations; nsteps=nsteps)
    
    # Plot results
    plot_results(truncations, cpu_times, gpu_times)
end

# Run the benchmarks
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
