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
    A = rand(spectral_grid.Grid, spectral_grid.nlat_half)   # same grid 
    B = rand(OctaHEALPixGrid, spectral_grid.nlat_half)      # different grid 
    f(lon, lat, sig) = sind(lon)*cosd(lat)*(1 - sig)

    # set things ...
    set!(simulation, vor=L, lf = lf)
    @test simulation.prognostic_variables.vor[lf] == L

    set!(simulation, div=L, lf = lf; add=true)
    @test simulation.prognostic_variables.div[lf] == (simulation.prognostic_variables.div[lf] .+ L)

    set!(simulation, temp=L2, lf=lf)

    set!(simulation, humid=L3, lf=lf)

    set!(simulation, pres=L[:,1], lf=lf)

    


#     sph_data = [rand(LowerTriangularMatrix{spectral_grid.NF}, lmax+1, mmax+1) for i=1:nlev]

#     SpeedyWeather.set_vorticity!(P, sph_data)
#     SpeedyWeather.set_divergence!(P, sph_data)
#     SpeedyWeather.set_temperature!(P, sph_data)
#     SpeedyWeather.set_humidity!(P, sph_data)
#     SpeedyWeather.set_pressure!(P, sph_data[1])

#     for i=1:nlev
#         @test P.layers[i].timesteps[lf].vor == sph_data[i]
#         @test P.layers[i].timesteps[lf].div == sph_data[i]
#         @test P.layers[i].timesteps[lf].temp == sph_data[i]
#         @test P.layers[i].timesteps[lf].humid == sph_data[i]
#     end 
#     @test P.surface.timesteps[lf].pres == sph_data[1]

#     @test SpeedyWeather.get_vorticity(P) == sph_data
#     @test SpeedyWeather.get_divergence(P) == sph_data
#     @test SpeedyWeather.get_temperature(P) == sph_data
#     @test SpeedyWeather.get_humidity(P) == sph_data
#     @test SpeedyWeather.get_pressure(P) == sph_data[1]

#     SpeedyWeather.set_vorticity!(P, 0)
#     for i=1:nlev 
#         @test all(P.layers[i].timesteps[lf].vor .== 0)
#     end

#     grid_data = [transform(sph_data[i], M.spectral_transform) for i in eachindex(sph_data)]

#     SpeedyWeather.set_vorticity!(P, grid_data)
#     SpeedyWeather.set_divergence!(P, grid_data)
#     SpeedyWeather.set_temperature!(P, grid_data)
#     SpeedyWeather.set_humidity!(P, grid_data)
#     SpeedyWeather.set_pressure!(P, grid_data[1])

#     for i=1:nlev 
#         @test P.layers[i].timesteps[lf].vor ≈   sph_data[i]
#         @test P.layers[i].timesteps[lf].div ≈   sph_data[i]
#         @test P.layers[i].timesteps[lf].temp ≈  sph_data[i]
#         @test P.layers[i].timesteps[lf].humid ≈ sph_data[i]
#     end 
#     @test P.surface.timesteps[lf].pres ≈ sph_data[1]

#     grid_data = [transform(sph_data[i], M.spectral_transform) for i in eachindex(sph_data)]

#     SpeedyWeather.set_vorticity!(P, grid_data, M)
#     SpeedyWeather.set_divergence!(P, grid_data, M)
#     SpeedyWeather.set_temperature!(P, grid_data, M)
#     SpeedyWeather.set_humidity!(P, grid_data, M)
#     SpeedyWeather.set_pressure!(P, grid_data[1], M)

#     for i=1:nlev 
#         @test P.layers[i].timesteps[lf].vor   ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].div   ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].temp  ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].humid ≈ sph_data[i]
#     end 
#     @test P.surface.timesteps[lf].pres ≈ sph_data[1]

#     # test setting matrices 
#     spectral_grid = SpectralGrid(Grid=FullGaussianGrid)
#     initial_conditions = StartFromRest()
#     M = PrimitiveWetModel(; spectral_grid, initial_conditions)
#     simulation = initialize!(M)
#     P = simulation.prognostic_variables
 
#     nlev = M.geometry.nlev
#     lmax = M.spectral_transform.lmax
#     mmax = M.spectral_transform.mmax
#     lf = 1

#     grid_data = [transform(sph_data[i], M.spectral_transform) for i in eachindex(sph_data)]
#     matrix_data = [Matrix(grid_data[i]) for i in eachindex(grid_data)]

#     SpeedyWeather.set_vorticity!(P, matrix_data)
#     SpeedyWeather.set_divergence!(P, matrix_data)
#     SpeedyWeather.set_temperature!(P, matrix_data)
#     SpeedyWeather.set_humidity!(P, matrix_data)
#     SpeedyWeather.set_pressure!(P, matrix_data[1])

#     for i=1:nlev 
#         @test P.layers[i].timesteps[lf].vor   ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].div   ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].temp  ≈ sph_data[i]
#         @test P.layers[i].timesteps[lf].humid ≈ sph_data[i]
#     end 
#     @test P.surface.timesteps[lf].pres ≈ sph_data[1]

# end