@testset "GPU interpolation" begin
    arch = SpeedyWeather.GPU()
    grid_in = HEALPixGrid(40, arch)
    grid_out = FullClenshawGrid(30, arch)
    interp = RingGrids.interpolator(grid_out, grid_in)
    
    field_in = on_architecture(arch, rand(grid_in))
    field_out = on_architecture(arch, zeros(grid_out))
    RingGrids.interpolate!(field_out, field_in, interp)
    
    cpu_arch = SpeedyWeather.CPU()
    grid_in_cpu = on_architecture(cpu_arch, grid_in)
    grid_out_cpu = on_architecture(cpu_arch, grid_out)
    field_in_cpu = on_architecture(cpu_arch, field_in)
    field_out_cpu = on_architecture(cpu_arch, field_out)
    interp_cpu = RingGrids.interpolator(grid_out_cpu, grid_in_cpu)
    RingGrids.interpolate!(field_out_cpu, field_in_cpu, interp_cpu)
    @test on_architecture(cpu_arch, field_out) ≈ field_out_cpu
end


