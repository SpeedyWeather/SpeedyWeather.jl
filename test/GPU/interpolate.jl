@testset "GPU interpolation" begin
    arch = SpeedyWeather.GPU()
    grid_in = FullClenshawGrid(50, arch)
    grid_out = HEALPixGrid(26, arch)
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

    
    field_out_cpu = rand(grid_out_cpu)
    field_out_gpu = on_architecture(arch, field_out_cpu)
   
    SpeedyWeather.RingGrids.update_locator!(interp_cpu, field_out_cpu)
    SpeedyWeather.RingGrids.update_locator!(interp, field_out_gpu)

    @test interp_cpu.locator.js == on_architecture(cpu_arch, interp.locator.js)
    @test interp_cpu.locator.ij_as == on_architecture(cpu_arch, interp.locator.ij_as)
    @test interp_cpu.locator.ij_bs == on_architecture(cpu_arch, interp.locator.ij_bs)
    @test interp_cpu.locator.ij_cs == on_architecture(cpu_arch, interp.locator.ij_cs)
    @test interp_cpu.locator.ij_ds == on_architecture(cpu_arch, interp.locator.ij_ds)
    @test interp_cpu.locator.Δys == on_architecture(cpu_arch, interp.locator.Δys)
    @test interp_cpu.locator.Δabs == on_architecture(cpu_arch, interp.locator.Δabs)
    @test interp_cpu.locator.Δcds == on_architecture(cpu_arch, interp.locator.Δcds)
end


