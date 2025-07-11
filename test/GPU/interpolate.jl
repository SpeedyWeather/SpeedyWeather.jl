@testset "GPU interpolation" begin
    arch = SpeedyWeather.GPU()
    grid_in = HEALPixGrid(40, architecture=arch)
    grid_out = FullClenshawGrid(30, architecture=arch)
    interp = RingGrids.interpolator(grid_out, grid_in)
    
    field_in = on_architecture(arch, rand(grid_in))
    field_out = on_architecture(arch, zeros(grid_out))
    interpolate!(field_out, field_in, interp)
    
    cpu_arch = SpeedyWeather.CPU()
    field_in_cpu = on_architecture(cpu_arch, field_in)
    field_out_cpu = on_architecture(cpu_arch, field_out)
    field_out_cpu = interpolate(field_out_cpu, field_in_cpu)
    @test on_architecture(cpu_arch, field_out) ≈ field_out_cpu
end