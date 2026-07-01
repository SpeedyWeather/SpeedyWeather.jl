@testset "Test Variables set!" begin

    nlayers = 8
    nlayers_soil = 2
    trunc = 31
    NF = Float64
    complex_NF = Complex{NF}
    spectral_grid = SpectralGrid(; NF, trunc, nlayers)  # define resolution
    model = PrimitiveWetModel(spectral_grid)            # construct model
    simulation = initialize!(model)                     # initialize all model components

    step = 2

    # test data
    L = rand(spectral_grid.SpectralVariable3D, trunc + 2, trunc + 1, nlayers)
    L_grid = transform(L, model.spectral_transform)

    L2 = rand(spectral_grid.SpectralVariable3D, trunc - 5, trunc - 6, nlayers)      # smaller
    L2_trunc = SpeedyTransforms.spectral_truncation(L2, size(L, 1, ZeroBased, as = Matrix), size(L, 2, ZeroBased, as = Matrix))
    L3 = rand(spectral_grid.SpectralVariable3D, trunc + 6, trunc + 5, nlayers)      # bigger
    L3_trunc = SpeedyTransforms.spectral_truncation(L3, size(L, 1, ZeroBased, as = Matrix), size(L, 2, ZeroBased, as = Matrix))

    A = rand(NF, spectral_grid.grid, nlayers)                                   # same grid
    A_spec = transform(A, model.spectral_transform)
    B = rand(NF, OctaHEALPixGrid, spectral_grid.nlat_half, nlayers)             # different grid
    D = rand(spectral_grid.GridVariable3D, spectral_grid.grid, nlayers_soil)    # 3D land data

    f(lon, lat, sig) = sind(lon) * cosd(lat) * (1 - sig)

    prog_new = simulation.variables.prognostic
    prog_old = deepcopy(prog_new)

    # set things ...

    # LTA
    set!(simulation, vorticity = L; step)
    @test get_step(prog_new.vorticity, step) == L

    set!(simulation, divergence = L; step, add = true)
    @test get_step(prog_new.divergence, step) == (get_step(prog_old.divergence, step) .+ L)

    set!(simulation, temperature = L2; step)
    @test get_step(prog_new.temperature, step) ≈ L2_trunc

    set!(simulation, humidity = L3; step)
    @test get_step(prog_new.humidity, step) ≈ L3_trunc

    set!(simulation, pressure = L[:, 1]; step)
    @test get_step(prog_new.pressure, step) == L[:, 1]

    set!(simulation, vorticity = A; step)
    @test get_step(prog_new.vorticity, step) == A_spec

    set!(simulation, vorticity = A; step, add = true)
    @test get_step(prog_new.vorticity, step) == (2 .* A_spec)

    # grids
    set!(simulation, sea_surface_temperature = A[:, 1], namespace = :ocean)
    @test get_step(prog_new.ocean.sea_surface_temperature, 1) == A[:, 1]

    set!(simulation, sea_ice_concentration = B[:, 1], namespace = :ocean, add = true)
    set!(simulation, sea_ice_concentration = B[:, 1], namespace = :ocean, add = true)
    C = similar(A[:, 1])
    RingGrids.interpolate!(C, B[:, 1]; NF)
    @test all(isapprox(prog_new.ocean.sea_ice_concentration, prog_old.ocean.sea_ice_concentration .+ C, atol = 1.0e-6))

    Di = deepcopy(prog_new.land.soil_temperature)
    RingGrids.interpolate!(Di, D; NF)
    set!(simulation, soil_temperature = D, namespace = :land)
    @test prog_new.land.soil_temperature == Di

    set!(simulation, soil_moisture = D; step, namespace = :land)
    @test prog_new.land.soil_moisture == Di

    # numbers
    set!(simulation, vorticity = Float32(3.0); step)
    M3 = zeros(NF, spectral_grid.grid, nlayers) .+ 3    # same grid
    M3_spec = transform(M3, model.spectral_transform)
    @test get_step(prog_new.vorticity, step) ≈ M3_spec

    set!(simulation, vorticity = Float32(3.0); step, add = true)
    @test get_step(prog_new.vorticity, step) ≈ (2 .* M3_spec)

    set!(simulation, sea_surface_temperature = Float16(3.0), namespace = :ocean)
    @test all(prog_new.ocean.sea_surface_temperature .≈ 3.0)

    set!(simulation, sea_surface_temperature = Float16(3.0), add = true, namespace = :ocean)
    @test all(prog_new.ocean.sea_surface_temperature .≈ 6.0)

    # vor_div, create u,v first in spectral space
    u = randn(spectral_grid.SpectralVariable3D, trunc + 2, trunc + 1, nlayers)
    v = randn(spectral_grid.SpectralVariable3D, trunc + 2, trunc + 1, nlayers)

    # set imaginary component of m=0 to 0 as the rotation of zonal modes is arbitrary
    SpeedyTransforms.zero_imaginary_zonal_modes!(u)
    SpeedyTransforms.zero_imaginary_zonal_modes!(v)

    SpeedyTransforms.spectral_truncation!(u, 25)     # truncate to lowest 11 wavenumbers
    SpeedyTransforms.spectral_truncation!(v, 25)

    u_grid = transform(u, model.spectral_transform)
    v_grid = transform(v, model.spectral_transform)

    set!(simulation, u = u_grid, v = v_grid, coslat_scaling_included = false; step)

    # now obtain U, V (scaled with coslat) from vorticity, div
    U = similar(u)
    V = similar(v)

    SpeedyTransforms.UV_from_vordiv!(U, V, get_step(prog_new.vorticity, step), get_step(prog_new.divergence, step), model.spectral_transform)

    # back to grid and unscale on the fly
    u_grid2 = transform(U, model.spectral_transform, unscale_coslat = true)
    v_grid2 = transform(V, model.spectral_transform, unscale_coslat = true)

    u2 = transform(u_grid2, model.spectral_transform)
    v2 = transform(v_grid2, model.spectral_transform)

    # for lm in eachindex(u, v, u2, v2)
    #     @test u[lm] ≈ u2[lm] atol = sqrt(sqrt(eps(spectral_grid.NF)))
    #     @test v[lm] ≈ v2[lm] atol = sqrt(sqrt(eps(spectral_grid.NF)))
    # end

    # functions
    (; londs, latds, σ_levels_full) = model.geometry
    for k in RingGrids.eachlayer(A)
        for ij in RingGrids.eachgridpoint(A)
            A[ij, k] = f(londs[ij], latds[ij], σ_levels_full[k])
        end
    end
    transform!(A_spec, A, model.spectral_transform)

    set!(simulation, vorticity = f; step)
    @test get_step(prog_new.vorticity, step) ≈ A_spec

    # groups
    set!(simulation, geopotential = 1, group = :dynamics)
    @test all(simulation.variables.dynamics.geopotential .== 1)
end

@testset "Set! grids" begin
    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            FullClenshawGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
        )
        spectral_grid = SpectralGrid(; Grid)
        (; GridVariable2D, nlat_half) = spectral_grid

        grid = zeros(GridVariable2D, nlat_half)

        @test all(set!(grid, 3) .== 3)
        set!(grid, (λ, φ) -> cosd(φ))
        @test all(0 .<= grid .<= 1)

        # with geometry
        geometry = Geometry(spectral_grid)
        set!(grid, (λ, φ) -> cosd(φ), geometry)
        @test all(0 .<= grid .<= 1)
    end
end

@testset "Set! albedo" begin
    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            FullClenshawGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
        )

        spectral_grid = SpectralGrid(; Grid)
        geometry = Geometry(spectral_grid)
        albedo = ManualAlbedo(spectral_grid)

        @test all(set!(albedo, 0.0625) .== 0.0625)
        set!(albedo.albedo, (λ, φ) -> clamp(sind(λ) * cosd(φ), 0.1, 1))
        @test all(0.1 .<= albedo.albedo .<= 1)

        # with geometry
        set!(albedo, (λ, φ) -> cosd(φ), geometry)
        @test all(0 .<= albedo.albedo .<= 1)
    end
end
