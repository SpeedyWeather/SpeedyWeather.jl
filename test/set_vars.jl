@testset "Test PrognosticVariables set_vars! and get_var" begin 

    # test setting LowerTriangularMatrices
    P, D, M = initialize_speedy(initial_conditions=:rest,model=PrimitiveEquation)

    nlev = M.geometry.nlev
    lmax = M.spectral_transform.lmax
    mmax = M.spectral_transform.mmax
    lf = 1

    sph_data = [rand(LowerTriangularMatrix, lmax+2, mmax+1) for i=1:nlev]

    SpeedyWeather.set_vorticity!(P, sph_data)
    SpeedyWeather.set_divergence!(P, sph_data)
    SpeedyWeather.set_temperature!(P, sph_data)
    SpeedyWeather.set_humidity!(P, sph_data)
    SpeedyWeather.set_pressure!(P, sph_data[1])

    for i=1:nlev 
        @test P.layers[i].leapfrog[lf].vor == sph_data[i]
        @test P.layers[i].leapfrog[lf].div == sph_data[i]
        @test P.layers[i].leapfrog[lf].temp == sph_data[i]
        @test P.layers[i].leapfrog[lf].humid == sph_data[i]
    end 
    @test P.pres.leapfrog[lf] == sph_data[1]

    @test SpeedyWeather.get_vorticity(P) == sph_data
    @test SpeedyWeather.get_divergence(P) == sph_data
    @test SpeedyWeather.get_temperature(P) == sph_data
    @test SpeedyWeather.get_humidity(P) == sph_data
    @test SpeedyWeather.get_pressure(P) == sph_data[1]

    SpeedyWeather.set_vorticity!(P, 0)
    for i=1:nlev 
        @test all(P.layers[i].leapfrog[lf].vor .== 0)
    end

    # test setting grids 
    P, D, M = initialize_speedy(initial_conditions=:rest,model=PrimitiveEquation)

    grid_data = [gridded(sph_data[i], M.spectral_transform) for i in eachindex(sph_data)]

    SpeedyWeather.set_vorticity!(P, grid_data)
    SpeedyWeather.set_divergence!(P, grid_data)
    SpeedyWeather.set_temperature!(P, grid_data)
    SpeedyWeather.set_humidity!(P, grid_data)
    SpeedyWeather.set_pressure!(P, grid_data[1])

    for i=1:nlev 
        @test all(isapprox(P.layers[i].leapfrog[lf].vor, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].div, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].temp, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].humid, sph_data[i]))
    end 
    @test all(isapprox(P.pres.leapfrog[lf],sph_data[1]))

    P, D, M = initialize_speedy(initial_conditions=:rest,model=PrimitiveEquation)

    grid_data = [gridded(sph_data[i], M.spectral_transform) for i in eachindex(sph_data)]

    SpeedyWeather.set_vorticity!(P, grid_data, M)
    SpeedyWeather.set_divergence!(P, grid_data, M)
    SpeedyWeather.set_temperature!(P, grid_data, M)
    SpeedyWeather.set_humidity!(P, grid_data, M)
    SpeedyWeather.set_pressure!(P, grid_data[1], M)

    for i=1:nlev 
        @test all(isapprox(P.layers[i].leapfrog[lf].vor, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].div, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].temp, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].humid, sph_data[i]))
    end 
    @test all(isapprox(P.pres.leapfrog[lf],sph_data[1]))

    # test setting matrices 
    P, D, M = initialize_speedy(initial_conditions=:rest,model=PrimitiveEquation)

    matrix_data = [Matrix(grid_data[i]) for i in eachindex(grid_data)]

    SpeedyWeather.set_vorticity!(P, grid_data)
    SpeedyWeather.set_divergence!(P, grid_data)
    SpeedyWeather.set_temperature!(P, grid_data)
    SpeedyWeather.set_humidity!(P, grid_data)
    SpeedyWeather.set_pressure!(P, grid_data[1])

    for i=1:nlev 
        @test all(isapprox(P.layers[i].leapfrog[lf].vor, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].div, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].temp, sph_data[i]))
        @test all(isapprox(P.layers[i].leapfrog[lf].humid, sph_data[i]))
    end 
    @test all(isapprox(P.pres.leapfrog[lf],sph_data[1]))

end 