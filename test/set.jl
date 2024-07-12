
import SpeedyWeather: set! 

@testset "Test PrognosticVariables set!" begin 

    N_lev = 8 
    N_trunc = 31
    spectral_grid = SpectralGrid(trunc=N_trunc, nlev=N_lev)          # define resolution
    model = PrimitiveWetModel(; spectral_grid)   # construct model
    simulation = initialize!(model)                         # initialize all model components
 
    lmax = model.spectral_transform.lmax
    mmax = model.spectral_transform.mmax
    lf = 2

    # test data 
    L = rand(LowerTriangularArray{Float32}, N_trunc+2, N_trunc+1, N_lev)
    L_grid = transform(L, model.spectral_transform)
    
    L2 = rand(LowerTriangularArray{Float32}, N_trunc-5, N_trunc-6, N_lev)  # smaller  
    L2_trunc = spectral_truncation(L2, size(L, 1, as=Matrix), size(L, 2, as=Matrix))
    
    L3 = rand(LowerTriangularArray{Float32}, N_trunc+6, N_trunc+5, N_lev)  # bigger 
    L3_trunc = spectral_truncation(L3, size(L, 1, as=Matrix), size(L, 2, as=Matrix))
    
    A = rand(spectral_grid.Grid{Float32}, spectral_grid.nlat_half, N_lev)   # same grid 
    A_spec = transform(A, model.spectral_transform)
    B = rand(OctaHEALPixGrid{Float32}, spectral_grid.nlat_half, N_lev)      # different grid 
    
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

    set!(simulation, land_surface_temperature=L[:,1], lf=lf, add=true)
    @test prog_new.land.land_surface_temperature == (prog_old.land.land_surface_temperature .+ L_grid[:,1])

    set!(simulation, snow_depth=L2[:,1], lf=lf)
    @test prog_new.land.snow_depth == L2_trunc[:,1]

    set!(simulation, soil_moisture_layer1=L3[:,1], lf=lf)
    @test prog_new.land.soil_moisture_layer1 == L3_trunc[:,1]

    set!(simulation, soil_moisture_layer2=L[:,1], lf=lf)   
    @test prog_new.land.soil_moisture_layer2 == L[:,1]

    # numbers
    set!(simulation, vor=Float16(3.), lf=lf)
    M3 = zeros(spectral_grid.Grid{Float32}, spectral_grid.nlat_half, N_lev) .+ 3  # same grid 
    M3_spec = transform(M3, model.spectral_transform)
    @test prog_new.vor[lf] ≈ M3_spec

    set!(simulation, vor=Float16(3.), lf=lf, add=true)
    @test prog_new.vor[lf] ≈ (2 .* M3_spec)

    set!(simulation, sea_surface_temperature=Float16(3.))
    @test all(prog_new.sea_surface_temperature .≈ 3.)

    set!(simulation, sea_surface_temperature=Float16(3.), add=true)
    @test all(prog_new.sea_surface_temperature .≈ 6.)

    # vor_div 
    set_vordiv!(simulation, u=A , v=(2 .* A))

    u2 = similar(A)
    v2 = similar(B)
    UV_from_vordiv!(u2, v2, prog_new.vor[lf], prog_new.div[lf],  )

    @test A ≈ u2 
    @test (2 .* A) ≈ v2

    # functions 
    (; londs, latds, σ_levels_full) = model.geometry
    for k in eachgrid(A)
        for ij in eachgridpoint(A)
            A[ij,k] = f(londs[ij],lats[ij],σ_levels_full[k])
        end 
    end 
    transform!(A_spec, A, model.spectral_transform)

    set!(simulation, vor=f; lf)
    @test prog_new.vor[lf] ≈ A_spec

end