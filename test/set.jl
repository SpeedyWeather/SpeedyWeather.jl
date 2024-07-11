@testset "Test PrognosticVariables set!" begin 

    N_lev = 8 
    N_trunc = 31
    spectral_grid = SpectralGrid(trunc=N_trunc, nlev=N_lev)          # define resolution
    model = PrimitiveWetModel(; spectral_grid)   # construct model
    initial_conditions = 
    simulation = initialize!(model)                         # initialize all model components
 
    lmax = M.spectral_transform.lmax
    mmax = M.spectral_transform.mmax
    lf = 2

    # test data 
    L = rand(LowerTriangularArray, N_trunc+2, N_trunc+1, N_lev)
    L2 = rand(LowerTriangularArray, N_trunc-5, N_trunc-6, N_lev)  # smaller  
    L3 = rand(LowerTriangularArray, N_trunc+6, N_trunc+5, N_lev)  # bigger 
    A = rand(spectral_grid.Grid, spectral_grid.nlat_half, N_lev)   # same grid 
    A_spec = transform(A, model.spectral_transform)
    B = rand(OctaHEALPixGrid, spectral_grid.nlat_half, N_lev)      # different grid 
    f(lon, lat, sig) = sind(lon)*cosd(lat)*(1 - sig)

    # set things ...

    # LTA
    set!(simulation, vor=L, lf = lf)
    @test simulation.prognostic_variables.vor[lf] == L

    set!(simulation, div=L, lf = lf; add=true)
    @test simulation.prognostic_variables.div[lf] == (simulation.prognostic_variables.div[lf] .+ L)

    set!(simulation, temp=L2, lf=lf)
    @test simulation.prognostic_variables.temp[lf] == 

    set!(simulation, humid=L3, lf=lf)
    
    set!(simulation, pres=L[:,1], lf=lf)
    @test simulation.prognostic_variables.pres[lf] == L[:,1]

    set!(simulation, vor=A, lf=lf)
    @test simulation.prognostic_variables.vor[lf] == A_spec

    set!(simulation, vor=A, lf=lf, add=true)
    @test simulation.prognostic_variables.vor[lf] == (simulation.prognostic_variables.vor[lf] .+ A_spec)

    # grids 
    set!(simulation, sea_surface_temperature=A[:,1], lf=lf)
    @test simulation.ocean.sea_surface_temperature == A[:,1]

    set!(simulation, sea_ice_concentration=B[:,1], lf=lf, add=true)



    set!(simulation, land_surface_temperature=L[:,1], lf=lf, add=true)

    set!(simulation, snow_depth=L2[:,1], lf=lf)

    set!(simulation, soil_moisture_layer1=L3[:,1], lf=lf)

    set!(simulation, soil_moisture_layer2=L[:,1], lf=lf)   

    # numbers
    set!(simulation, vor=Float16(3.), lf=lf)
    @test all(simulation.vor[lf])

    set!(simulation, vor=Float16(3.), lf=lf, add=true)

    set!(simulation, sea_surface_temperature=Float16(3.))

    set!(simulation, sea_surface_temperature=Float16(3.), add=true)

    # functions 



    # vor_div 


# end