"""
Test script to validate GPU kernel implementations for implicit corrections.

This script verifies that the GPU kernel versions produce the same results
as the CPU versions within numerical tolerance.
"""

using SpeedyWeather
using Test, BenchmarkTools

@testset "Implicit Kernel Tests" begin
    
    @testset "Shallow Water Implicit Correction" begin
        println("\nTesting Shallow Water Implicit Correction...")
        
        NF = Float32
        trunc = 31
        nlayers = 1
        
        # Setup CPU model
        spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, architecture=CPU())
        model_cpu = ShallowWaterModel(; spectral_grid=spectral_grid_cpu)
        sim_cpu = initialize!(model_cpu)
        
        # Run one timestep to get realistic data
        run!(sim_cpu, steps=20)
        
        # Copy state for GPU test
        diagn_cpu = sim_cpu.diagnostic_variables
        progn_cpu = sim_cpu.prognostic_variables
        
        # Save original tendencies
        div_tend_orig = copy(diagn_cpu.tendencies.div_tend)
        pres_tend_orig = copy(diagn_cpu.tendencies.pres_tend)
        
        # Run CPU version
        SpeedyWeather.implicit_correction_cpu!(diagn_cpu, progn_cpu, model_cpu.implicit, model_cpu)
        
        div_tend_cpu = copy(diagn_cpu.tendencies.div_tend)
        pres_tend_cpu = copy(diagn_cpu.tendencies.pres_tend)
        
        # Reset tendencies
        diagn_cpu.tendencies.div_tend .= div_tend_orig
        diagn_cpu.tendencies.pres_tend .= pres_tend_orig
        
        # Run GPU version (on CPU architecture, kernel still works)
        SpeedyWeather.implicit_correction!(diagn_cpu, progn_cpu, model_cpu.implicit, model_cpu)
        
        div_tend_gpu = diagn_cpu.tendencies.div_tend
        pres_tend_gpu = diagn_cpu.tendencies.pres_tend
        
        # Compare results
        @test div_tend_cpu ≈ div_tend_gpu rtol=1e-5
        @test pres_tend_cpu ≈ pres_tend_gpu rtol=1e-5
        
        println("  ✓ Shallow Water kernels produce correct results")
    end
    
    @testset "Primitive Equation Implicit Correction" begin
        println("\nTesting Primitive Equation Implicit Correction...")
        
        NF = Float32
        trunc = 63
        nlayers = 8
        
        # Setup CPU model
        spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, architecture=CPU())
        model_cpu = PrimitiveDryModel(; spectral_grid=spectral_grid_cpu)
        sim_cpu = initialize!(model_cpu)
        
        run!(sim_cpu, steps=10)
        diagn_cpu = sim_cpu.diagnostic_variables
        progn_cpu = sim_cpu.prognostic_variables

        fill!(diagn_cpu.tendencies, 0, PrimitiveDry)
        SpeedyWeather.dynamics_tendencies!(diagn_cpu, progn_cpu, 2, model_cpu)

        # ===== TEST 1: CPU VERSION =====
        # Make complete deep copies of ALL inputs before running
        diagn_cpu_copy = deepcopy(diagn_cpu)
        progn_cpu_copy = deepcopy(progn_cpu)
        model_cpu_copy = deepcopy(model_cpu)
        
        #@btime SpeedyWeather.implicit_correction_cpu!($diagn_cpu_copy, $progn_cpu_copy, $model_cpu_copy.implicit, $model_cpu_copy)
        SpeedyWeather.implicit_correction_cpu!(diagn_cpu_copy, progn_cpu_copy, model_cpu_copy.implicit, model_cpu_copy)

        div_tend_cpu = copy(diagn_cpu_copy.tendencies.div_tend)
        temp_tend_cpu = copy(diagn_cpu_copy.tendencies.temp_tend)
        pres_tend_cpu = copy(diagn_cpu_copy.tendencies.pres_tend)
        
        # ===== TEST 2: 4-KERNEL GPU VERSION =====
        # Fresh deep copies for GPU version
        diagn_gpu_copy = deepcopy(diagn_cpu)
        progn_gpu_copy = deepcopy(progn_cpu)
        model_gpu_copy = deepcopy(model_cpu)
        
        #@btime SpeedyWeather.implicit_correction!($diagn_gpu_copy, $progn_gpu_copy, $model_gpu_copy.implicit, $model_gpu_copy)
        SpeedyWeather.implicit_correction_kernels!(diagn_gpu_copy, progn_gpu_copy, model_gpu_copy.implicit, model_gpu_copy)

        div_tend_gpu = copy(diagn_gpu_copy.tendencies.div_tend)
        temp_tend_gpu = copy(diagn_gpu_copy.tendencies.temp_tend)
        pres_tend_gpu = copy(diagn_gpu_copy.tendencies.pres_tend)
        
        # Compare results
        @test div_tend_cpu ≈ div_tend_gpu rtol=1e-5
        @test temp_tend_cpu ≈ temp_tend_gpu rtol=1e-5
        @test pres_tend_cpu ≈ pres_tend_gpu rtol=1e-5
        
        println("  ✓ Primitive Equation kernels (4-kernel version) produce correct results")
        
        # ===== TEST 3: SINGLE-KERNEL VERSION =====
        # Fresh deep copies for single-kernel version
        diagn_lm_copy = deepcopy(diagn_cpu)
        progn_lm_copy = deepcopy(progn_cpu)
        model_lm_copy = deepcopy(model_cpu)
        
        #@btime SpeedyWeather.implicit_correction_lm_only!($diagn_lm_copy, $progn_lm_copy, $model_lm_copy.implicit, $model_lm_copy)
        SpeedyWeather.implicit_correction_lm_only!(diagn_lm_copy, progn_lm_copy, model_lm_copy.implicit, model_lm_copy)

        div_tend_lm = diagn_lm_copy.tendencies.div_tend
        temp_tend_lm = diagn_lm_copy.tendencies.temp_tend
        pres_tend_lm = diagn_lm_copy.tendencies.pres_tend
        
        # Compare single-kernel results with CPU
        @test div_tend_cpu ≈ div_tend_lm rtol=1e-5
        @test temp_tend_cpu ≈ temp_tend_lm rtol=1e-5
        @test pres_tend_cpu ≈ pres_tend_lm rtol=1e-5
        
        println("  ✓ Primitive Equation kernels (single-kernel version) produce correct results")
    end
    
    # Test with actual GPU if available
    if GPU in SpeedyWeather.Architectures.available_architectures()
        @testset "GPU vs CPU Comparison" begin
            println("\nTesting GPU vs CPU comparison...")
            
            NF = Float32
            trunc = 31
            nlayers = 8
            
            # Setup CPU model
            spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, architecture=CPU())
            model_cpu = PrimitiveDryModel(; spectral_grid=spectral_grid_cpu)
            sim_cpu = initialize!(model_cpu)
            timestep!(sim_cpu, 1)
            
            # Setup GPU model
            spectral_grid_gpu = SpectralGrid(; NF, trunc, nlayers, architecture=GPU())
            model_gpu = PrimitiveDryModel(; spectral_grid=spectral_grid_gpu)
            sim_gpu = initialize!(model_gpu)
            timestep!(sim_gpu, 1)
            
            # Transfer CPU data to GPU
            diagn_gpu = sim_gpu.diagnostic_variables
            progn_gpu = sim_gpu.prognostic_variables
            
            # Copy tendencies from CPU to GPU
            diagn_gpu.tendencies.div_tend .= on_architecture(GPU(), sim_cpu.diagnostic_variables.tendencies.div_tend)
            diagn_gpu.tendencies.temp_tend .= on_architecture(GPU(), sim_cpu.diagnostic_variables.tendencies.temp_tend)
            diagn_gpu.tendencies.pres_tend .= on_architecture(GPU(), sim_cpu.diagnostic_variables.tendencies.pres_tend)
            
            # Run CPU version
            SpeedyWeather.implicit_correction_cpu!(
                sim_cpu.diagnostic_variables, 
                sim_cpu.prognostic_variables, 
                model_cpu.implicit, 
                model_cpu
            )
            
            # Run GPU version
            SpeedyWeather.implicit_correction!(diagn_gpu, progn_gpu, model_gpu.implicit, model_gpu)
            synchronize(model_gpu.spectral_grid.architecture)
            
            # Transfer GPU results back to CPU for comparison
            div_tend_cpu = sim_cpu.diagnostic_variables.tendencies.div_tend
            div_tend_gpu = on_architecture(CPU(), diagn_gpu.tendencies.div_tend)
            
            temp_tend_cpu = sim_cpu.diagnostic_variables.tendencies.temp_tend
            temp_tend_gpu = on_architecture(CPU(), diagn_gpu.tendencies.temp_tend)
            
            pres_tend_cpu = sim_cpu.diagnostic_variables.tendencies.pres_tend
            pres_tend_gpu = on_architecture(CPU(), diagn_gpu.tendencies.pres_tend)
            
            # Compare results
            @test div_tend_cpu ≈ div_tend_gpu rtol=1e-4
            @test temp_tend_cpu ≈ temp_tend_gpu rtol=1e-4
            @test pres_tend_cpu ≈ pres_tend_gpu rtol=1e-4
            
            println("  ✓ GPU and CPU produce consistent results")
        end
    else
        println("\nGPU not available, skipping GPU-specific tests")
    end
end

println("\n" * "="^70)
println("All implicit kernel tests passed!")
println("="^70)
