using CUDA 
using SpeedyWeather
using BenchmarkTools 
using CairoMakie

# Define the range of array sizes (N) to benchmark
array_sizes = [31, 63, 127, 255, 511]
nlayers_list = [1, 8, 32, 64]
float_types = [Float32]

# Updated function to generate random inputs with device parameter
function generate_random_inputs(N, nlayers, T, device)
    spectral_grid = SpectralGrid(NF=T, trunc=N, nlayers=nlayers, device=device)
    S = SpectralTransform(spectral_grid)
    grids = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
    specs = rand(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
    if device == SpeedyWeather.GPU()
        specs = cu(specs)
        grids = cu(grids)
    end
    return S, specs, grids
end

# Single run_benchmarks function with device parameter
function run_benchmarks(array_sizes, nlayers_list, float_types, device)
    results = Dict{String, Dict{DataType, Matrix{Float64}}}()

    for T in float_types
        # Initialize times matrices for each function
        # times_legendre_specs_S = Matrix{Float64}(undef, length(array_sizes), length(nlayers_list))
        times_legendre_backward = Matrix{Float64}(undef, length(array_sizes), length(nlayers_list))
        times_fourier_forward = Matrix{Float64}(undef, length(array_sizes), length(nlayers_list))
        times_fourier_backward = Matrix{Float64}(undef, length(array_sizes), length(nlayers_list))

        for (j, nlayers) in enumerate(nlayers_list)
            for (i, N) in enumerate(array_sizes)
                # Generate inputs
                S, specs, grids = generate_random_inputs(N, nlayers, T, device)
                
                println("Running benchmark for device=$device, N=$N, nlayers=$nlayers, T=$T")
                println()

                # # Time forward legendre
                # time = @elapsed SpeedyTransforms._legendre!(specs, S.scratch_memory_north, S.scratch_memory_south, S)
                # times_legendre_specs_S[i, j] = time

                # Time inverse legendre
                time = @elapsed SpeedyTransforms._legendre!(S.scratch_memory_north, S.scratch_memory_south, specs, S)
                times_legendre_backward[i, j] = time

                # Time _fourier!(..., grids, S)
                time = @elapsed SpeedyTransforms._fourier!(S.scratch_memory_north, S.scratch_memory_south, grids, S)
                times_fourier_backward[i, j] = time
                
                # Time _fourier!(grids, ..., S)
                time = @elapsed SpeedyTransforms._fourier!(grids, S.scratch_memory_north, S.scratch_memory_south, S)
                times_fourier_forward[i, j] = time
            end
        end

        # Store results for each function
        # results["forward_legendre"] = get(results, "forward_legendre", Dict{DataType, Matrix{Float64}}())
        # results["forward_legendre"][T] = times_legendre_specs_S

        results["inverse_legendre"] = get(results, "inverse_legendre", Dict{DataType, Matrix{Float64}}())
        results["inverse_legendre"][T] = times_legendre_backward

        results["forward_fourier"] = get(results, "forward_fourier", Dict{DataType, Matrix{Float64}}())
        results["forward_fourier"][T] = times_fourier_forward

        results["inverse_fourier"] = get(results, "inverse_fourier", Dict{DataType, Matrix{Float64}}())
        results["inverse_fourier"][T] = times_fourier_backward
    end
    return results
end

function plot_speedup(cpu_results, gpu_results)
    # Create a larger figure with subplots
    fig = Figure(size = (800, 1600))
    n = length(nlayers_list)
    
    for (j, nlayers) in enumerate(nlayers_list)
        # Create a subplot for each nlayers value
        ax = Axis(fig[j, 1], xlabel="Array Size (N)", ylabel="Speedup (CPU/GPU)", title="nlayers = $nlayers", width=400, height=300, yscale=log10)

        for (func_name, type_results) in cpu_results
            for (T, cpu_times) in type_results
                gpu_times = gpu_results[func_name][T]
                speedup = cpu_times ./ gpu_times
                # Convert array_sizes and speedup to a vector of Point2 objects
                points = [Point2(array_sizes[i], speedup[i, j]) for i in 1:length(array_sizes)]
                lines!(ax, points, label="$func_name ($T)", linewidth=2)
            end
        end

        # Add a horizontal line for y = 1
        max_size = maximum(array_sizes)
        lines!(ax, [0, max_size], [1, 1], linestyle=:dot, color=:black, label="y = 1")
        axislegend(ax)
    end

    # Save the figure
    save("benchmark_speedup.png", fig)
end

function plot_times(cpu_results, gpu_results)
    # Create a larger figure with subplots
    fig = Figure(size = (800, 1600))
    n = length(nlayers_list)
    
    for (j, nlayers) in enumerate(nlayers_list)
        # Create a subplot for each nlayers value
        ax = Axis(fig[j, 1], xlabel="Array Size (N)", ylabel="Speedup (CPU/GPU)", title="nlayers = $nlayers", width=400, height=300)

        for (results_vec, label) in zip([cpu_results, gpu_results], ["CPU", "GPU"])
            for (func_name, type_results) in results_vec
                for (T, times) in type_results
                    # Convert array_sizes and times to a vector of Point2 objects
                    points = [Point2(array_sizes[i], times[i, j]) for i in 1:length(array_sizes)]
                    lines!(ax, points, label="$func_name ($label)", linewidth=2)
                end
            end
        end

       axislegend(ax)
    end

    # Save the figure
    save("benchmark_times.png", fig)
end

# Run benchmarks for CPU and GPU
@show cpu_results = run_benchmarks(array_sizes, nlayers_list, float_types, SpeedyWeather.CPU())
@show gpu_results = run_benchmarks(array_sizes, nlayers_list, float_types, SpeedyWeather.GPU())

# Plot the results
plot_speedup(cpu_results, gpu_results)
plot_times(cpu_results, gpu_results)
