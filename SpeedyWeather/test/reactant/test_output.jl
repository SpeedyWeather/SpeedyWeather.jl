# Output writers can't interpolate on Reactant arrays outside of compiled code, so they
# interpolate on the CPU for a Reactant model, moving the model's fields there first.

@testset "Output writers interpolate on the CPU for Reactant" begin
    arch = ReactantDevice()
    cpu = SpeedyWeather.CPU()
    @test SpeedyWeather.output_architecture(arch) isa SpeedyWeather.CPU

    spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION)
    spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION)
    output = NetCDFOutput(spectral_grid)
    output_cpu = NetCDFOutput(spectral_grid_cpu)

    @test SpeedyWeather.ismatching(cpu, output.field2D.data)
    @test SpeedyWeather.ismatching(cpu, output.interpolator.locator.ij_as)
    @test output.host2D === output.field2D          # written from directly, no extra copy

    # a field on the Reactant device interpolates to the same as on the CPU
    field_cpu = rand(Float32, spectral_grid_cpu.grid)
    field = on_architecture(arch, field_cpu)
    SpeedyWeather.interpolate_output!(output, output.field2D, field)
    SpeedyWeather.interpolate_output!(output_cpu, output_cpu.field2D, field_cpu)
    @test output.field2D.data == output_cpu.field2D.data
end
