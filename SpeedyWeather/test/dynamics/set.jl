@testset "Test Variables set!" begin

    nlayers = 8
    nlayers_soil = 2
    truncation = 32
    NF = Float64
    complex_NF = Complex{NF}
    spectral_grid = SpectralGrid(; NF, truncation, nlayers)     # define resolution
    model = PrimitiveWetModel(spectral_grid)                    # construct model
    simulation = initialize!(model)                             # initialize all model components

    step = 2
    # ocean and land are stepped with EulerForward by default (1 step), use their last (=current) step
    ocean_step = SpeedyWeather.nsteps(simulation.variables.prognostic.ocean.sea_surface_temperature)
    land_step = SpeedyWeather.nsteps(simulation.variables.prognostic.land.soil_temperature)

    # test data
    L = rand(spectral_grid.SpectralVariable3D, truncation + 1, truncation, nlayers)
    L_grid = transform(L, model.spectral_transform)

    L2 = rand(spectral_grid.SpectralVariable3D, truncation - 5, truncation - 6, nlayers)      # smaller
    L2_trunc = SpeedyTransforms.spectral_truncation(L2, size(L, 1, as = Matrix), size(L, 2, as = Matrix))
    L3 = rand(spectral_grid.SpectralVariable3D, truncation + 6, truncation + 5, nlayers)      # bigger
    L3_trunc = SpeedyTransforms.spectral_truncation(L3, size(L, 1, as = Matrix), size(L, 2, as = Matrix))

    # set! removes the extra degree lmax+1 that scalar (=all prognostic) variables don't use
    L_trunc = SpeedyTransforms.spectral_truncation!(copy(L))
    SpeedyTransforms.spectral_truncation!(L3_trunc)

    A = rand(NF, spectral_grid.grid, nlayers)                                   # same grid
    A_spec = transform(A, model.spectral_transform)
    A_spec_trunc = SpeedyTransforms.spectral_truncation!(copy(A_spec))
    B = rand(NF, OctaHEALPixGrid, spectral_grid.nlat_half, nlayers)             # different grid
    D = rand(spectral_grid.GridVariable3D, spectral_grid.grid, nlayers_soil)    # 3D land data

    f(lon, lat, sig) = sind(lon) * cosd(lat) * (1 - sig)

    prog_new = simulation.variables.prognostic
    prog_old = deepcopy(prog_new)

    # set things ...

    # LTA
    set!(simulation, vorticity = L; step, truncate = false)     # keep the extra degree
    @test get_step(prog_new.vorticity, step) == L

    set!(simulation, vorticity = L; step)
    @test get_step(prog_new.vorticity, step) == L_trunc

    set!(simulation, divergence = L; step, add = true)
    @test get_step(prog_new.divergence, step) == (get_step(prog_old.divergence, step) .+ L_trunc)

    set!(simulation, temperature = L2; step)
    @test get_step(prog_new.temperature, step) ≈ L2_trunc

    set!(simulation, humidity = L3; step)
    @test get_step(prog_new.humidity, step) ≈ L3_trunc

    set!(simulation, pressure = L[:, 1]; step)
    @test get_step(prog_new.pressure, step) == L_trunc[:, 1]

    # set! transforms in-place into the variable (a view), A_spec into a new array: the `muladd`s in
    # the Legendre transform may be fused differently for these, ≈ instead of == (differs on x86)
    set!(simulation, vorticity = A; step)
    @test get_step(prog_new.vorticity, step) ≈ A_spec_trunc

    set!(simulation, vorticity = A; step, add = true)
    @test get_step(prog_new.vorticity, step) ≈ (2 .* A_spec_trunc)

    # grids
    set!(simulation, sea_surface_temperature = A[:, 1]; namespace = :ocean)
    @test prog_new.ocean.sea_surface_temperature[:, 1] == A[:, 1]
    set!(simulation, sea_surface_temperature = A[:, 1]; step = ocean_step, namespace = :ocean)
    @test prog_new.ocean.sea_surface_temperature[:, ocean_step] == A[:, 1]

    set!(simulation, sea_ice_concentration = B[:, 1]; step = ocean_step, namespace = :ocean, add = true)
    C = similar(A[:, 1])
    RingGrids.interpolate!(C, B[:, 1]; NF)

    sic_new = get_step(prog_new.ocean.sea_ice_concentration, ocean_step)
    sic_old = get_step(prog_old.ocean.sea_ice_concentration, ocean_step)
    @test all(isapprox(sic_new, sic_old .+ C, atol = 1.0e-6))

    lst = get_step(prog_new.land.soil_temperature, land_step)
    Di = deepcopy(lst)
    RingGrids.interpolate!(Di, D; NF)
    set!(simulation, soil_temperature = D, namespace = :land, step = land_step)
    @test lst == Di

    set!(simulation, soil_moisture = D; namespace = :land)
    @test get_step(prog_new.land.soil_moisture, 1) == Di

    # numbers
    set!(simulation, vorticity = Float32(3.0); step)
    M3 = zeros(NF, spectral_grid.grid, nlayers) .+ 3    # same grid
    M3_spec = transform(M3, model.spectral_transform)
    @test get_step(prog_new.vorticity, step) ≈ M3_spec

    set!(simulation, vorticity = Float32(3.0); step, add = true)
    @test get_step(prog_new.vorticity, step) ≈ (2 .* M3_spec)

    set!(simulation, sea_surface_temperature = Float16(3.0), namespace = :ocean)
    @test all(prog_new.ocean.sea_surface_temperature[:, 1] .≈ 3.0)

    set!(simulation, sea_surface_temperature = Float16(3.0), step = 1, add = true, namespace = :ocean)
    @test all(prog_new.ocean.sea_surface_temperature[:, 1] .≈ 6.0)

    set!(simulation, sea_surface_temperature = 5, step = ocean_step, namespace = :ocean)
    @test all(prog_new.ocean.sea_surface_temperature .== 5)

    # vor_div, create u,v first in spectral space
    u = randn(spectral_grid.SpectralVariable3D, truncation + 1, truncation, nlayers)
    v = randn(spectral_grid.SpectralVariable3D, truncation + 1, truncation, nlayers)

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
    SpeedyTransforms.spectral_truncation!(A_spec)

    set!(simulation, vorticity = f; step)
    @test get_step(prog_new.vorticity, step) ≈ A_spec

    # groups
    # non-stepped variables need step = nothing
    set!(simulation, geopotential = 1, group = :dynamics)
    @test all(simulation.variables.dynamics.geopotential .== 1)
end

@testset "Set! u, v" begin
    nlayers = 8
    truncation = 32
    NF = Float64
    spectral_grid = SpectralGrid(; NF, truncation, nlayers)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    S = model.spectral_transform
    (; londs, latds, σ_levels_full) = model.geometry
    radius = model.planet.radius

    vor = get_step(simulation.variables.prognostic.vorticity, 1)
    div = get_step(simulation.variables.prognostic.divergence, 1)

    # solid body rotation u = U₀cos(φ), v = 0 has vorticity 2U₀sin(φ)/radius and no divergence
    U₀ = 20
    set!(simulation, u = (λ, φ, σ) -> U₀ * cosd(φ), v = (λ, φ, σ) -> 0)
    vor_grid = transform(vor, S)
    for k in 1:nlayers
        @test vor_grid.data[:, k] ≈ 2U₀ * sind.(latds) / radius
    end
    @test maximum(abs, div.data) < sqrt(eps(NF)) * U₀ / radius

    # random large-scale u, v in spectral space and on the grid
    u = randn(spectral_grid.SpectralVariable3D, truncation + 1, truncation, nlayers)
    v = randn(spectral_grid.SpectralVariable3D, truncation + 1, truncation, nlayers)
    SpeedyTransforms.zero_imaginary_zonal_modes!(u)
    SpeedyTransforms.zero_imaginary_zonal_modes!(v)
    SpeedyTransforms.spectral_truncation!(u, 25)
    SpeedyTransforms.spectral_truncation!(v, 25)
    u_grid = transform(u, S)
    v_grid = transform(v, S)

    set!(simulation, u = u_grid, v = v_grid)
    vor_grid_path, div_grid_path = copy(vor), copy(div)
    @test !all(iszero, vor_grid_path.data)
    @test !all(iszero, div_grid_path.data)

    # same from spectral u, v (v was previously computed from u)
    set!(simulation, u = u, v = v)
    @test vor ≈ vor_grid_path
    @test div ≈ div_grid_path

    # same without spectral transform provided
    set!(simulation.variables, model.geometry; u = u_grid, v = v_grid)
    @test vor ≈ vor_grid_path
    @test div ≈ div_grid_path

    # functions are evaluated on the grid, no spectral round trip before the 1/coslat scaling
    f_u(λ, φ, σ) = 10cosd(φ) + 5sind(λ) * cosd(φ)^2 * σ
    f_v(λ, φ, σ) = 3cosd(λ) * cosd(φ) * sind(φ)
    u_grid_f = zeros(NF, spectral_grid.grid, nlayers)
    v_grid_f = zeros(NF, spectral_grid.grid, nlayers)
    for k in RingGrids.eachlayer(u_grid_f), ij in RingGrids.eachgridpoint(u_grid_f)
        u_grid_f[ij, k] = f_u(londs[ij], latds[ij], σ_levels_full[k])
        v_grid_f[ij, k] = f_v(londs[ij], latds[ij], σ_levels_full[k])
    end

    set!(simulation, u = u_grid_f, v = v_grid_f)
    vor_f, div_f = copy(vor), copy(div)
    set!(simulation, u = f_u, v = f_v)
    @test vor ≈ vor_f
    @test div ≈ div_f

    # smaller spectral u, v are zero-padded (the size check was a chained comparison)
    u_small = SpeedyTransforms.spectral_truncation(u, 21, 20)
    v_small = SpeedyTransforms.spectral_truncation(v, 21, 20)
    u_pad = SpeedyTransforms.spectral_truncation(u_small, truncation + 1, truncation)
    v_pad = SpeedyTransforms.spectral_truncation(v_small, truncation + 1, truncation)
    for coslat_scaling_included in (false, true)
        set!(simulation; u = u_pad, v = v_pad, coslat_scaling_included)
        vor_pad, div_pad = copy(vor), copy(div)
        set!(simulation; u = u_small, v = v_small, coslat_scaling_included)
        @test vor ≈ vor_pad
        @test div ≈ div_pad
    end

    # add
    set!(simulation, u = u_grid, v = v_grid)
    set!(simulation, u = u_grid, v = v_grid, add = true)
    @test vor ≈ 2 .* vor_grid_path
    @test div ≈ 2 .* div_grid_path

    set!(simulation, u = 0, v = 0, add = true)      # adding nothing
    @test vor ≈ 2 .* vor_grid_path
    @test div ≈ 2 .* div_grid_path

    set!(simulation, u = 0, v = 0)
    @test all(iszero, vor.data)
    @test all(iszero, div.data)
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
