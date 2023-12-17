@testset "Geopotential reasonable" begin
    for NF in (Float32,Float64)
        nlev = 8
        spectral_grid = SpectralGrid(;NF,nlev,Grid=FullGaussianGrid)
        model = PrimitiveWetModel(;spectral_grid)
        simulation = initialize!(model)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        # give every layer some constant temperature
        temp = 280      # in Kelvin
        lf = 1
        for (progn_layer,diagn_layer) in zip(p.layers,d.layers)
            progn_layer_lf = progn_layer.timesteps[lf]
            progn_layer_lf.temp[1] = temp*model.spectral_transform.norm_sphere
            fill!(progn_layer_lf.humid,0)                               # dry core
            SpeedyWeather.gridded!(diagn_layer,progn_layer_lf,model)    # propagate spectral state to grid
            SpeedyWeather.linear_virtual_temperature!(diagn_layer,progn_layer,model,lf)
        end

        SpeedyWeather.geopotential!(d,model.orography,model.constants)

        # approximate heights [m] for this setup
        heights = [27000,18000,13000,9000,6000,3700,1800,700] 

        for k in 1:8
            geopot_grid = Matrix(gridded(d.layers[k].dynamics_variables.geopot))
            height_over_ocean = geopot_grid[48,24]/model.planet.gravity    # middle of pacific
            @test heights[k] ≈ height_over_ocean rtol=0.5                         # very large error allowed
        end
    end
end

@testset "Add geopotential and kinetic energy, compute -∇²B term, no errors" begin
    for NF in (Float32,Float64)
        nlev = 8
        spectral_grid = SpectralGrid(;NF,nlev,Grid=FullGaussianGrid)
        m = PrimitiveWetModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        # give every layer some constant temperature
        temp = 280      # in Kelvin
        for k in 1:nlev
            p.layers[k].timesteps[1].temp[1] = temp*m.spectral_transform.norm_sphere
        end

        SpeedyWeather.geopotential!(d,m.orography,m.constants)

        lf = 1
        for (progn_layer,diagn_layer) in zip(p.layers,d.layers)
            progn_layer_lf = progn_layer.timesteps[lf]
            SpeedyWeather.gridded!(diagn_layer,progn_layer_lf,m)             # propagate spectral state to grid
            SpeedyWeather.bernoulli_potential!(diagn_layer,m.spectral_transform)
        end
    end
end

@testset "Virtual temperature calculation" begin
    for NF in (Float32,Float64)
        nlev = 8
        spectral_grid = SpectralGrid(;NF,nlev,Grid=FullGaussianGrid)
        m = PrimitiveWetModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        # give every layer some constant temperature
        temp = 280      # in Kelvin
        lf = 1
        for (progn_layer,diagn_layer) in zip(p.layers,d.layers)
            progn_layer_lf = progn_layer.timesteps[lf]
            progn_layer_lf.temp[1] = temp*m.spectral_transform.norm_sphere
            fill!(progn_layer_lf.humid,0)                 
            progn_layer_lf.humid[1] = rand(NF)
            SpeedyWeather.gridded!(diagn_layer,progn_layer_lf,m)    # propagate spectral state to grid
            
            temp_lf = progn_layer.timesteps[lf].temp     # always passed on for compat with DryCore
            SpeedyWeather.virtual_temperature!(diagn_layer,temp_lf,m)

            T  = diagn_layer.grid_variables.temp_grid
            Tv = diagn_layer.grid_variables.temp_virt_grid

            for ij in SpeedyWeather.eachgridpoint(T,Tv)
                @test T[ij] < Tv[ij]    # virtual temperature has to be higher as q > 0
            end
        end
    end
end