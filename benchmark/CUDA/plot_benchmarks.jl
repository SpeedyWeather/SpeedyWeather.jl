using CUDA 
using SpeedyWeather
using Adapt
using BenchmarkTools 
using CairoMakie

include("benchmark_suite.jl")

# Define the range of array sizes (N) to benchmark
trunc_list = [15, 31, 63, 127, 255, 511]
nlayers_list = [1, 8, 32, 64]
float_types = [Float32]
grid_list = [FullGaussianGrid]

# # Updated function to generate random inputs with device parameter
# function generate_random_inputs(N, nlayers, T, device)
#     spectral_grid = SpectralGrid(NF=T, trunc=N, nlayers=nlayers, device=device)
#     S = SpectralTransform(spectral_grid)
#     grids = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
#     specs = rand(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
#     if device == SpeedyWeather.GPU()
#         specs = cu(specs)
#         grids = cu(grids)
#     end
#     return S, specs, grids
# end

function get_median_result(trial::BenchmarkTools.Trial)
    t = median(trial)

    return t.time
end


# Single run_benchmarks function with device parameter
function run_benchmarks(trunc_list, nlayers_list, float_types, device)
    results = Dict{String, Dict{DataType, Matrix{Float64}}}()
    Grid = grid_list[1]

    for NF in float_types
        # Initialize times matrices for each function
        times_legendre_forward = Matrix{Float64}(undef, length(trunc_list), length(nlayers_list))
        times_legendre_backward = Matrix{Float64}(undef, length(trunc_list), length(nlayers_list))
        times_fourier_forward = Matrix{Float64}(undef, length(trunc_list), length(nlayers_list))
        times_fourier_backward = Matrix{Float64}(undef, length(trunc_list), length(nlayers_list))

        for (j, nlayers) in enumerate(nlayers_list)
            for (i, trunc) in enumerate(trunc_list)
                # Generate inputs
                spectral_grid = SpectralGrid(;NF, trunc, Grid, nlayers, device)
                S = SpectralTransform(spectral_grid)
                specs, grids = generate_random_inputs(spectral_grid)
                # S, specs, grids = generate_random_inputs(N, nlayers, T, device)
                
                println("Running benchmark for device=$device, trunc=$trunc, nlayers=$nlayers, NF=$NF")
                println()

                # Time forward legendre
                b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($specs, $S.scratch_memory_north, $S.scratch_memory_south, $S)
                times_legendre_forward[i, j] = median(b).time

                # Time inverse legendre
                b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($S.scratch_memory_north, $S.scratch_memory_south, $specs, $S)
                times_legendre_backward[i, j] = median(b).time

                # Time _fourier!(..., grids, S)
                b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($S.scratch_memory_north, $S.scratch_memory_south, $grids, $S)
                times_fourier_backward[i, j] = median(b).time
                
                # Time _fourier!(grids, ..., S)
                b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($grids, $S.scratch_memory_north, $S.scratch_memory_south, $S)
                times_fourier_forward[i, j] = median(b).time
            end
        end

        # Store results for each function
        results["forward_legendre"] = get(results, "forward_legendre", Dict{DataType, Matrix{Float64}}())
        results["forward_legendre"][NF] = times_legendre_forward

        results["inverse_legendre"] = get(results, "inverse_legendre", Dict{DataType, Matrix{Float64}}())
        results["inverse_legendre"][NF] = times_legendre_backward

        results["forward_fourier"] = get(results, "forward_fourier", Dict{DataType, Matrix{Float64}}())
        results["forward_fourier"][NF] = times_fourier_forward

        results["inverse_fourier"] = get(results, "inverse_fourier", Dict{DataType, Matrix{Float64}}())
        results["inverse_fourier"][NF] = times_fourier_backward
    end
    return results
end

function plot_speedup(cpu_results, gpu_results, figx=500, figy=1000)
    # Create a larger figure with subplots
    fig = Figure(size = (figy, figx))
    n = length(nlayers_list)
    axes = Vector{Axis}(undef, n)

    # Create a subplot for each nlayers value
    for (j, nlayers) in enumerate(nlayers_list)
        axes[j] = Axis(fig[1, j], xlabel="Trunc (T)", title="nlayers = $nlayers", yscale=log10) 
        ax = axes[j]
        
        for (func_name, type_results) in cpu_results
            for (T, cpu_times) in type_results
                gpu_times = gpu_results[func_name][T]
                speedup = cpu_times ./ gpu_times
                # Convert array_sizes and speedup to a vector of Point2 objects
                points = [Point2(trunc_list[i], speedup[i, j]) for i in 1:length(trunc_list)]
                lines!(ax, points, label="$func_name", linewidth=2)
            end
        end

        # Add a horizontal line for y = 1
        max_size = maximum(trunc_list)
        lines!(ax, [0, max_size], [1, 1], linestyle=:dot, color=:black, label="y = 1")
                
        if j == n
            fig[1, n+1] = Legend(fig, ax, "Function", framevisible=false)            
        end
    end

    axes[1].ylabel = "Speedup (CPU/GPU)"
    for ax in axes[2:n]
        linkyaxes!(axes[1], ax)
        ax.yticklabelsvisible = false 
    end

    # Save the figure
    save("benchmark_speedup.png", fig)
end

function plot_times(cpu_results, gpu_results, figx=500, figy=1000)
    # Create a larger figure with subplots
    fig = Figure(size = (figy, figx))
    n = length(nlayers_list)
    axes = Vector{Axis}(undef, n)
    
    # Create a subplot for each nlayers value
    for (j, nlayers) in enumerate(nlayers_list)
        axes[j] = Axis(fig[1, j], xlabel="Trunc (T)", title="nlayers = $nlayers", yscale=log10) 
        ax = axes[j]
        
        for (i, (func_name, type_results)) in enumerate(cpu_results)    
            cpu_times = cpu_results[func_name][Float32] ./ 1e9
            gpu_times = gpu_results[func_name][Float32] ./ 1e9
            cpu_points = [Point2(trunc_list[i], cpu_times[i, j]) for i in 1:length(trunc_list)]
            gpu_points = [Point2(trunc_list[i], gpu_times[i, j]) for i in 1:length(trunc_list)]
            
            # Convert array_sizes and speedup to a vector of Point2 objects
            lines!(ax, cpu_points, label="$func_name (CPU)", linewidth=2, 
                linestyle=:solid, color=i, colormap=:tab10, colorrange=(1,10))
            lines!(ax, gpu_points, label="$func_name (GPU)", linewidth=2, 
                linestyle=:dash, color=i, colormap=:tab10, colorrange=(1,10))
        end

        if j == n
            fig[1, n+1] = Legend(fig, ax, "Function", framevisible=false)            
        end
    end

    axes[1].ylabel = "Time elapsed (s)"
    for ax in axes[2:n]
        linkyaxes!(axes[1], ax)
        ax.yticklabelsvisible = false 
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
