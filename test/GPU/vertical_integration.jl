# The vertical_integration! function is
# one of the few functions that have distinct GPU version. 
# Here, we just test that they are identical to the CPU version
# The CPU version is in turn tested more properly in the regular unit test

@testset "Vertical Integration" begin
    arch_gpu = SpeedyWeather.GPU()
    arch_cpu = SpeedyWeather.CPU()
    spectral_grid_cpu = SpectralGrid(architecture=arch_cpu)
    spectral_grid_gpu = SpectralGrid(architecture=arch_gpu)

    # don't initialize parameterizations
    model_cpu = PrimitiveWetModel(spectral_grid_cpu, physics=false)
    model_gpu = PrimitiveWetModel(spectral_grid_gpu, physics=false)

    simulation_cpu = initialize!(model_cpu)
    simulation_gpu = initialize!(model_gpu)

    run!(simulation_cpu, steps=2)

    progn_cpu, diagn_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
    progn_gpu, diagn_gpu, model_gpu = SpeedyWeather.unpack(simulation_gpu)

    # copy all values needed to have have them equal
    diagn_gpu.dynamics.∇lnp_x .= on_architecture(arch_gpu, diagn_cpu.dynamics.∇lnp_x)
    diagn_gpu.dynamics.∇lnp_y .= on_architecture(arch_gpu, diagn_cpu.dynamics.∇lnp_y)

    diagn_gpu.grid.u_grid = on_architecture(arch_gpu, diagn_cpu.grid.u_grid)
    diagn_gpu.grid.v_grid = on_architecture(arch_gpu, diagn_cpu.grid.v_grid)
    diagn_gpu.grid.div_grid = on_architecture(arch_gpu, diagn_cpu.grid.div_grid)

    diagn_gpu.dynamics.u_mean_grid = on_architecture(arch_gpu, diagn_cpu.dynamics.u_mean_grid)
    diagn_gpu.dynamics.v_mean_grid = on_architecture(arch_gpu, diagn_cpu.dynamics.v_mean_grid)
    diagn_gpu.dynamics.div_mean_grid = on_architecture(arch_gpu, diagn_cpu.dynamics.div_mean_grid)
    diagn_gpu.dynamics.div_mean = on_architecture(arch_gpu, diagn_cpu.dynamics.div_mean)

    diagn_gpu.dynamics.div_sum_above = on_architecture(arch_gpu, diagn_cpu.dynamics.div_sum_above)
    diagn_gpu.dynamics.uv∇lnp_sum_above = on_architecture(arch_gpu, diagn_cpu.dynamics.uv∇lnp_sum_above)

    progn_gpu.div .= on_architecture(arch_gpu, progn_cpu.div)

    # Run vertical integration on both CPU and GPU
    SpeedyWeather.vertical_integration!(diagn_cpu, progn_cpu, 1, model_cpu.geometry)
    SpeedyWeather.vertical_integration!(diagn_gpu, progn_gpu, 1, model_gpu.geometry)

    # Check that results are identical
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.u_mean_grid) ≈ diagn_cpu.dynamics.u_mean_grid
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.v_mean_grid) ≈ diagn_cpu.dynamics.v_mean_grid
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.div_mean_grid) ≈ diagn_cpu.dynamics.div_mean_grid
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.div_mean) ≈ diagn_cpu.dynamics.div_mean
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.div_sum_above) ≈ diagn_cpu.dynamics.div_sum_above
    @test on_architecture(arch_cpu, diagn_gpu.dynamics.uv∇lnp_sum_above) ≈ diagn_cpu.dynamics.uv∇lnp_sum_above
end