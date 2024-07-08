@testset "Geopotential reasonable" begin
    for NF in (Float32, Float64)
        nlev = 8
        spectral_grid = SpectralGrid(; NF, nlev, Grid=FullGaussianGrid)
        model = PrimitiveWetModel(; spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables
        temp = progn.temp[1]
        humid = progn.humid[1]

        # give every layer some constant temperature by setting the l=m=0 mode (index 1) for all k
        temp0 = 280      # in Kelvin
        temp[1, :] .= temp0 * model.spectral_transform.norm_sphere
        humid .= 0

        lf = 1      # leapfrog time step
        SpeedyWeather.transform!(diagn, progn, lf, model)
        SpeedyWeather.linear_virtual_temperature!(diagn, progn, lf, model)
        SpeedyWeather.geopotential!(diagn, model.geopotential, model.orography)
        geopot_grid = Array(transform(diagn.dynamics.geopot, model.spectral_transform))
        
        # approximate heights [m] for this setup
        heights = [27000, 18000, 13000, 9000, 6000, 3700, 1800, 700] 

        height_over_ocean = geopot_grid[48, 24, :]/model.planet.gravity    # middle of pacific
        for k in eachmatrix(temp)
            @test heights[k] ≈ height_over_ocean[k] rtol=0.5                         # very large error allowed
        end
    end
end

@testset "Add geopotential and kinetic energy, compute -∇²B term, no errors" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, nlayers=8, Grid=FullGaussianGrid)
        model = PrimitiveWetModel(; spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables
        temp = progn.temp[1]
        humid = progn.humid[1]

        # give every layer some constant temperature by setting the l=m=0 mode (index 1) for all k
        temp0 = 280      # in Kelvin
        temp[1, :] .= temp0 * model.spectral_transform.norm_sphere
        humid .= 0

        lf = 1      # leapfrog time step
        SpeedyWeather.transform!(diagn, progn, lf, model)
        SpeedyWeather.linear_virtual_temperature!(diagn, progn, lf, model)
        SpeedyWeather.geopotential!(diagn, model.geopotential, model.orography)
        SpeedyWeather.bernoulli_potential!(diagn, model.spectral_transform)
    end
end

@testset "Virtual temperature calculation" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, nlayers=8, Grid=FullGaussianGrid)
        model = PrimitiveWetModel(; spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables

        temp = progn.temp[1]
        humid = progn.humid[1]

        # give every layer some constant temperature by setting the l=m=0 mode (index 1) for all k
        temp0 = 280      # in Kelvin
        temp .= 0
        temp[1, :] .= temp0 * model.spectral_transform.norm_sphere
        
        humid .= 0
        humid[1, :] .= 1e-2*(rand(spectral_grid.nlayers) .+ 0.1) * model.spectral_transform.norm_sphere

        lf = 1      # leapfrog time step
        SpeedyWeather.transform!(diagn, progn, lf, model)
        SpeedyWeather.virtual_temperature!(diagn, model)

        T  = diagn.grid.temp_grid
        Tv = diagn.grid.temp_virt_grid

        for ijk in eachindex(T, Tv)
            # virtual temperature has to be higher as q > 0
            @test T[ijk] < Tv[ijk]
        end
    end
end