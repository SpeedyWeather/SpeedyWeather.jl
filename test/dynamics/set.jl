@testset "Test PrognosticVariables set!" begin 

    nlayers = 8
    nlayers_soil = 2
    trunc = 31
    NF = Float64
    complex_NF = Complex{NF}
    spectral_grid = SpectralGrid(; NF, trunc, nlayers)  # define resolution
    model = PrimitiveWetModel(spectral_grid)            # construct model
    simulation = initialize!(model)                     # initialize all model components
 
    lf = 2

    # test data 
    L = rand(spectral_grid.SpectralVariable3D, trunc+2, trunc+1, nlayers)
    L_grid = transform(L, model.spectral_transform)
    
    L2 = rand(spectral_grid.SpectralVariable3D, trunc-5, trunc-6, nlayers)      # smaller  
    L2_trunc = spectral_truncation(L2, size(L, 1, ZeroBased, as=Matrix), size(L, 2, ZeroBased, as=Matrix))
    L3 = rand(spectral_grid.SpectralVariable3D, trunc+6, trunc+5, nlayers)      # bigger 
    L3_trunc = spectral_truncation(L3, size(L, 1, ZeroBased, as=Matrix), size(L, 2, ZeroBased, as=Matrix))
    
    A = rand(NF, spectral_grid.grid, nlayers)                                   # same grid 
    A_spec = transform(A, model.spectral_transform)
    B = rand(NF, OctaHEALPixGrid, spectral_grid.nlat_half, nlayers)             # different grid 
    D = rand(spectral_grid.GridVariable3D, spectral_grid.grid, nlayers_soil)    # 3D land data

    f(lon, lat, sig) = sind(lon)*cosd(lat)*(1 - sig)

    prog_new = simulation.prognostic_variables
    prog_old = deepcopy(prog_new)

    # set things ...

    # LTA
    set!(simulation, vor=L, lf = lf)
    @test get_step(prog_new.vor, lf) == L

    set!(simulation, div=L, lf = lf; add=true)
    @test get_step(prog_new.div, lf) == (get_step(prog_old.div, lf) .+ L)

    set!(simulation, temp=L2, lf=lf)
    @test get_step(prog_new.temp, lf) ≈ L2_trunc

    set!(simulation, humid=L3, lf=lf)
    @test get_step(prog_new.humid, lf) ≈ L3_trunc
    
    set!(simulation, pres=L[:,1], lf=lf)
    @test get_step(prog_new.pres, lf) == L[:,1]

    set!(simulation, vor=A, lf=lf)
    @test get_step(prog_new.vor, lf) == A_spec

    set!(simulation, vor=A, lf=lf, add=true)
    @test get_step(prog_new.vor, lf) == (2 .* A_spec)

    # grids 
    set!(simulation, sea_surface_temperature=A[:,1], lf=lf)
    @test prog_new.ocean.sea_surface_temperature == A[:,1]

    set!(simulation, sea_ice_concentration=B[:,1], lf=lf, add=true)
    C = similar(A[:,1])
    SpeedyWeather.RingGrids.interpolate!(C, B[:,1])
    @test prog_new.ocean.sea_ice_concentration ≈ (prog_old.ocean.sea_ice_concentration .+ C)

    Di = deepcopy(prog_new.land.soil_temperature)
    RingGrids.interpolate!(Di, D)
    set!(simulation, soil_temperature=D, lf=lf)
    @test prog_new.land.soil_temperature == Di

    set!(simulation, soil_moisture=D, lf=lf)   
    @test prog_new.land.soil_moisture == Di

    # numbers
    set!(simulation, vor=Float32(3.), lf=lf)
    M3 = zeros(NF, spectral_grid.grid, nlayers) .+ 3    # same grid 
    M3_spec = transform(M3, model.spectral_transform)
    @test get_step(prog_new.vor, lf) ≈ M3_spec

    set!(simulation, vor=Float32(3.), lf=lf, add=true)
    @test get_step(prog_new.vor, lf) ≈ (2 .* M3_spec)

    set!(simulation, sea_surface_temperature=Float16(3.))
    @test all(prog_new.ocean.sea_surface_temperature .≈ 3.)

    set!(simulation, sea_surface_temperature=Float16(3.), add=true)
    @test all(prog_new.ocean.sea_surface_temperature .≈ 6.)

    # vor_div, create u,v first in spectral space
    u = randn(spectral_grid.SpectralVariable3D, trunc+2, trunc+1, nlayers)
    v = randn(spectral_grid.SpectralVariable3D, trunc+2, trunc+1, nlayers)
    
    # set imaginary component of m=0 to 0 as the rotation of zonal modes is arbitrary
    SpeedyTransforms.zero_imaginary_zonal_modes!(u)
    SpeedyTransforms.zero_imaginary_zonal_modes!(v)

    spectral_truncation!(u, 25)     # truncate to lowest 11 wavenumbers
    spectral_truncation!(v, 25)

    u_grid = transform(u, model.spectral_transform)
    v_grid = transform(v, model.spectral_transform)   

    set!(simulation, u=u_grid, v=v_grid, coslat_scaling_included=false, lf=lf)

    # now obtain U, V (scaled with coslat) from vor, div
    U = similar(u)
    V = similar(v)

    SpeedyTransforms.UV_from_vordiv!(U, V, get_step(prog_new.vor, lf), get_step(prog_new.div, lf), model.spectral_transform)

    # back to grid and unscale on the fly
    u_grid2 = transform(U, model.spectral_transform, unscale_coslat=true)
    v_grid2 = transform(V, model.spectral_transform, unscale_coslat=true)

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
            A[ij,k] = f(londs[ij], latds[ij], σ_levels_full[k])
        end 
    end 
    transform!(A_spec, A, model.spectral_transform)

    set!(simulation, vor=f; lf)
    @test get_step(prog_new.vor, lf) ≈ A_spec
end

@testset "Set! grids" begin
    @testset for Grid in (FullGaussianGrid,
                        OctahedralGaussianGrid,
                        FullClenshawGrid,
                        OctahedralClenshawGrid,
                        OctaminimalGaussianGrid)
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
    @testset for Grid in (FullGaussianGrid,
                        OctahedralGaussianGrid,
                        FullClenshawGrid,
                        OctahedralClenshawGrid,
                        OctaminimalGaussianGrid)

        spectral_grid = SpectralGrid(; Grid)
        geometry = Geometry(spectral_grid)
        albedo = ManualAlbedo(spectral_grid)

        @test all(set!(albedo, 0.0625) .== 0.0625)
        set!(albedo.albedo, (λ, φ) -> clamp(sind(λ)*cosd(φ), 0.1, 1))
        @test all(0.1 .<= albedo.albedo .<= 1)

        # with geometry
        set!(albedo, (λ, φ) -> cosd(φ), geometry)
        @test all(0 .<= albedo.albedo .<= 1)
    end
end