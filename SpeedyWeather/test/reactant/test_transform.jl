# The MatrixSpectralTransform is also called outside of compiled code, e.g. for the initial
# conditions, so it has to work eagerly on Reactant arrays in both directions

@testset "MatrixSpectralTransform on Reactant outside of compilation" begin
    arch = ReactantDevice()
    nlayers = 2
    spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION, nlayers)
    spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION, nlayers)
    M = MatrixSpectralTransform(spectral_grid)
    M_cpu = MatrixSpectralTransform(spectral_grid_cpu)

    # grid to spectral
    field_cpu = rand(Float32, spectral_grid_cpu.grid, nlayers)
    coeffs_cpu = transform(field_cpu, M_cpu)
    coeffs = transform(on_architecture(arch, field_cpu), M)
    @test on_architecture(SpeedyWeather.CPU(), coeffs.data) ≈ coeffs_cpu.data

    # spectral to grid
    field_back_cpu = transform(coeffs_cpu, M_cpu)
    field_back = transform(coeffs, M)
    @test on_architecture(SpeedyWeather.CPU(), field_back.data) ≈ field_back_cpu.data
end
