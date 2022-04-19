@testset "Spectral gradients no error test" begin
    for NF in (Float32,Float64)
        prog_vars,diag_vars,model_setup = SpeedyWeather.initialize_speedy(NF)

        temp_surf = prog_vars.temp[:,:,end]
        SpeedyWeather.gradient_longitude(temp_surf)
        SpeedyWeather.gradient_latitude(temp_surf)

        SpeedyWeather.∇²(temp_surf)
        SpeedyWeather.∇⁻²(temp_surf)

        temp_surf_grid = gridded(temp_surf)
        
        u_surf_grid = diag_vars.grid_variables.u_grid[:,:,end]
        v_surf_grid = diag_vars.grid_variables.v_grid[:,:,end]
        vor_surf_grid = diag_vars.grid_variables.vor_grid[:,:,end]

        SpeedyWeather.divergence_uvω_spectral(  u_surf_grid,
                                                v_surf_grid,
                                                vor_surf_grid,
                                                model_setup.geospectral)
    end
end