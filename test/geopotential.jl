@testset "Geopotential no errors" begin
    for NF in (Float32,Float64)
        nlev = 8
        p,d,m = initialize_speedy(NF,nlev=nlev,model=PrimitiveEquation,
                                        Grid=FullGaussianGrid)

        # give every layer some constant temperature
        temp = 280      # in Kelvin
        for k in 1:nlev
            p.layers[k].leapfrog[1].temp[1] = temp*m.spectral_transform.norm_sphere
        end

        SpeedyWeather.geopotential!(d,p,1,m.boundaries,m.geometry)

        # approximate heights [m] for this setup
        heights = [27000,18000,13000,9000,6000,3700,1800,700]

        for k in 1:8
            geopot_grid = Matrix(gridded(d.layers[k].dynamics_variables.geopot))
            height_over_ocean = geopot_grid[48,24]/m.parameters.gravity      # middle of pacific
            @test heights[k] ≈ height_over_ocean rtol=0.5                   # very large error allowed
        end
    end
end