# The vertical_integration! function is
# one of the few functions that have distinct GPU version.
# Here, we just test that they are identical to the CPU version
# The CPU version is in turn tested more properly in the regular unit test

@testset "Vertical Integration" begin
    arch_gpu = SpeedyWeather.GPU()
    arch_cpu = SpeedyWeather.CPU()

    spectral_grid_cpu = SpectralGrid(architecture = arch_cpu)
    spectral_grid_gpu = SpectralGrid(architecture = arch_gpu)

    # don't initialize parameterizations
    model_cpu = PrimitiveWetModel(spectral_grid_cpu, dynamics_only = true)
    model_gpu = PrimitiveWetModel(spectral_grid_gpu, dynamics_only = true)

    simulation_cpu = initialize!(model_cpu)
    simulation_gpu = initialize!(model_gpu)

    run!(simulation_cpu, steps = 2)

    vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
    vars_gpu, model_gpu = SpeedyWeather.unpack(simulation_gpu)

    # copy all values needed to have have them equal
    vars_gpu.dynamics.dpres_dx .= on_architecture(arch_gpu, vars_cpu.dynamics.dpres_dx)
    vars_gpu.dynamics.dpres_dy .= on_architecture(arch_gpu, vars_cpu.dynamics.dpres_dy)

    vars_gpu.grid.u .= on_architecture(arch_gpu, vars_cpu.grid.u)
    vars_gpu.grid.v .= on_architecture(arch_gpu, vars_cpu.grid.v)
    vars_gpu.grid.divergence .= on_architecture(arch_gpu, vars_cpu.grid.divergence)

    vars_gpu.dynamics.u_mean_grid .= on_architecture(arch_gpu, vars_cpu.dynamics.u_mean_grid)
    vars_gpu.dynamics.v_mean_grid .= on_architecture(arch_gpu, vars_cpu.dynamics.v_mean_grid)
    vars_gpu.dynamics.div_mean_grid .= on_architecture(arch_gpu, vars_cpu.dynamics.div_mean_grid)
    vars_gpu.dynamics.div_mean .= on_architecture(arch_gpu, vars_cpu.dynamics.div_mean)

    vars_gpu.dynamics.div_sum_above .= on_architecture(arch_gpu, vars_cpu.dynamics.div_sum_above)
    vars_gpu.dynamics.pres_flux_sum_above .= on_architecture(arch_gpu, vars_cpu.dynamics.pres_flux_sum_above)

    vars_gpu.prognostic.divergence .= on_architecture(arch_gpu, vars_cpu.prognostic.divergence)

    # Run vertical integration on both CPU and GPU
    SpeedyWeather.vertical_integration!(vars_cpu, 1, model_cpu.geometry)
    SpeedyWeather.vertical_integration!(vars_gpu, 1, model_gpu.geometry)

    # Check that results are identical
    @test on_architecture(arch_cpu, vars_gpu.dynamics.u_mean_grid) ≈ vars_cpu.dynamics.u_mean_grid
    @test on_architecture(arch_cpu, vars_gpu.dynamics.v_mean_grid) ≈ vars_cpu.dynamics.v_mean_grid
    @test on_architecture(arch_cpu, vars_gpu.dynamics.div_mean_grid) ≈ vars_cpu.dynamics.div_mean_grid
    @test on_architecture(arch_cpu, vars_gpu.dynamics.div_mean) ≈ vars_cpu.dynamics.div_mean
    @test on_architecture(arch_cpu, vars_gpu.dynamics.div_sum_above) ≈ vars_cpu.dynamics.div_sum_above
    @test on_architecture(arch_cpu, vars_gpu.dynamics.pres_flux_sum_above) ≈ vars_cpu.dynamics.pres_flux_sum_above
end
