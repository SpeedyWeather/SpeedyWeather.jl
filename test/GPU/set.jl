@testset "set! functions on GPU" begin 
    
    # 3D 
    architecture = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(; architecture)
    geometry = Geometry(spectral_grid)
    field = zeros(Float32, spectral_grid.grid, 8)

    set!(field, (λ, φ, η) -> 1f0 - η, geometry)

    @test field[1,:] ≈ (1f0 .- geometry.σ_levels_full)
    
    # zfac is just a pretty random function from Speedy with three input 
    # arguments, that's compiled and not dynamically generated
    set!(field, SpeedyWeather.zfac, geometry; enforce_static_func=true)
    
    # just compare to CPU version 
    architecture_cpu = SpeedyWeather.CPU()
    spectral_grid_cpu = SpectralGrid(; architecture_cpu)
    geometry_cpu = Geometry(spectral_grid_cpu)
    field_cpu = zeros(Float32, spectral_grid_cpu.grid, 8)

    set!(field_cpu, SpeedyWeather.zfac, geometry_cpu)

    @test field_cpu ≈ on_architecture(architecture_cpu, field)

    # 2D 
    field = zeros(Float32, spectral_grid.grid)
    set!(field, (λ, φ) -> φ, geometry)
    @test on_architecture(CPU(), field) ≈ RingGrids.get_londlatds(spectral_grid.grid)[2]

end