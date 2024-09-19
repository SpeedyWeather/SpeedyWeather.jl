@testset "Test PrognosticVariables set!" begin 

    nlayers = 8 
    trunc = 31
    spectral_grid = SpectralGrid(; trunc, nlayers)         # define resolution
    model = PrimitiveWetModel(; spectral_grid)              # construct model
    simulation = initialize!(model)                         # initialize all model components
 
    lmax = model.spectral_transform.lmax
    mmax = model.spectral_transform.mmax
    lf = 2

    # test data 
    L = rand(LowerTriangularArray{ComplexF32}, trunc+2, trunc+1, nlayers)
    L_grid = transform(L, model.spectral_transform)
    
    L2 = rand(LowerTriangularArray{ComplexF32}, trunc-5, trunc-6, nlayers)    # smaller  
    L2_trunc = spectral_truncation(L2, size(L, 1, ZeroBased, as=Matrix), size(L, 2, ZeroBased, as=Matrix))
    L3 = rand(LowerTriangularArray{ComplexF32}, trunc+6, trunc+5, nlayers)    # bigger 
    L3_trunc = spectral_truncation(L3, size(L, 1, ZeroBased, as=Matrix), size(L, 2, ZeroBased, as=Matrix))
    
    A = rand(spectral_grid.Grid{Float32}, spectral_grid.nlat_half, nlayers)       # same grid 
    A_spec = transform(A, model.spectral_transform)
    B = rand(OctaHEALPixGrid{Float32}, spectral_grid.nlat_half, nlayers)          # different grid 
    
    f(lon, lat, sig) = sind(lon)*cosd(lat)*(1 - sig)

    prog_new = simulation.prognostic_variables
    prog_old = deepcopy(prog_new)

    # set things ...

    # LTA
    set!(simulation, vor=L, lf = lf)
    @test prog_new.vor[lf] == L

    set!(simulation, div=L, lf = lf; add=true)
    @test prog_new.div[lf] == (prog_old.div[lf] .+ L)

    set!(simulation, temp=L2, lf=lf)
    @test prog_new.temp[lf] ≈ L2_trunc

    set!(simulation, humid=L3, lf=lf)
    @test prog_new.humid[lf] ≈ L3_trunc
    
    set!(simulation, pres=L[:,1], lf=lf)
    @test prog_new.pres[lf] == L[:,1]

    set!(simulation, vor=A, lf=lf)
    @test prog_new.vor[lf] == A_spec

    set!(simulation, vor=A, lf=lf, add=true)
    @test prog_new.vor[lf] == (2 .* A_spec)

    # grids 
    set!(simulation, sea_surface_temperature=A[:,1], lf=lf)
    @test prog_new.ocean.sea_surface_temperature == A[:,1]

    set!(simulation, sea_ice_concentration=B[:,1], lf=lf, add=true)
    C = similar(A[:,1])
    SpeedyWeather.RingGrids.interpolate!(C, B[:,1])
    @test prog_new.ocean.sea_ice_concentration ≈ (prog_old.ocean.sea_ice_concentration .+ C)

    set!(simulation, land_surface_temperature=L[:,1], lf=lf)
    @test prog_new.land.land_surface_temperature ≈ L_grid[:,1]

    set!(simulation, soil_moisture_layer2=L[:,1], lf=lf)   
    @test prog_new.land.soil_moisture_layer2 ≈ L_grid[:,1]

    # numbers
    set!(simulation, vor=Float32(3.), lf=lf)
    M3 = zeros(spectral_grid.Grid{Float32}, spectral_grid.nlat_half, nlayers) .+ 3  # same grid 
    M3_spec = transform(M3, model.spectral_transform)
    @test prog_new.vor[lf] ≈ M3_spec

    set!(simulation, vor=Float32(3.), lf=lf, add=true)
    @test prog_new.vor[lf] ≈ (2 .* M3_spec)

    set!(simulation, sea_surface_temperature=Float16(3.))
    @test all(prog_new.ocean.sea_surface_temperature .≈ 3.)

    set!(simulation, sea_surface_temperature=Float16(3.), add=true)
    @test all(prog_new.ocean.sea_surface_temperature .≈ 6.)

    # vor_div, create u,v first in spectral space
    u = randn(spectral_grid.SpectralVariable3D, trunc+2, trunc+1, nlayers)
    v = randn(spectral_grid.SpectralVariable3D, trunc+2, trunc+1, nlayers)
    
    u_grid = transform(u, model.spectral_transform)
    v_grid = transform(v, model.spectral_transform)   

    set!(simulation, u=u_grid, v=v_grid, coslat_scaling_included=false)

    # now obtain U, V (scaled with coslat) from vor, div
    U = similar(u)
    V = similar(v)

    SpeedyTransforms.UV_from_vordiv!(U, V, prog_new.vor[lf], prog_new.div[lf], model.spectral_transform)

    # back to grid and unscale on the fly
    u_grid2 = transform(U, model.spectral_transform, unscale_coslat=true)
    v_grid2 = transform(V, model.spectral_transform, unscale_coslat=true)

    # for ij in eachindex(u_grid, v_grid, u_grid2, v_grid2)
    #     @test_broken u_grid[ij] ≈ u_grid2[ij] 
    #     @test_broken v_grid[ij] ≈ v_grid2[ij]
    # end

    # functions 
    (; londs, latds, σ_levels_full) = model.geometry
    for k in SpeedyWeather.RingGrids.eachgrid(A)
        for ij in SpeedyWeather.RingGrids.eachgridpoint(A)
            A[ij,k] = f(londs[ij],latds[ij],σ_levels_full[k])
        end 
    end 
    transform!(A_spec, A, model.spectral_transform)

    set!(simulation, vor=f; lf)
    @test prog_new.vor[lf] ≈ A_spec
end