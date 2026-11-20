# a function of three arguments that is compiled and not dynamically generated, for `static_func = true`
static_test_function(λ, φ, η) = cosd(φ) * (1 - η) + sind(λ) / 10

@testset "set! functions on GPU" begin

    # 3D
    architecture = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(; architecture)
    geometry = Geometry(spectral_grid)
    field = zeros(Float32, spectral_grid.grid, 8)

    set!(field, (λ, φ, η) -> 1.0f0 - η, geometry, static_func = false)

    @test field[1, :] ≈ (1.0f0 .- geometry.σ_levels_full)

    set!(field, static_test_function, geometry; static_func = true)

    # just compare to CPU version
    architecture_cpu = SpeedyWeather.CPU()
    spectral_grid_cpu = SpectralGrid(; architecture = architecture_cpu)
    geometry_cpu = Geometry(spectral_grid_cpu)
    field_cpu = zeros(Float32, spectral_grid_cpu.grid, 8)

    set!(field_cpu, static_test_function, geometry_cpu)

    @test field_cpu ≈ on_architecture(architecture_cpu, field)

    # 2D
    field = zeros(Float32, spectral_grid.grid)
    set!(field, (λ, φ) -> φ, geometry)
    @test on_architecture(SpeedyWeather.CPU(), field) ≈ RingGrids.get_londlatds(spectral_grid.grid)[2]
end
