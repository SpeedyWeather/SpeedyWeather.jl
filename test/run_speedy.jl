@testset "Call run_speedy" begin
    run_speedy(Float64)     # just check that no error is triggered
    Prog = run_speedy(Float32)

    temp_surf = Prog.temp[:,:,end]
    SpeedyWeather.gradient_longitude(temp_surf)
    SpeedyWeather.gradient_latitude(temp_surf)

    SpeedyWeather.∇²(temp_surf)
    SpeedyWeather.∇⁻²(temp_surf)

    temp_surf_grid = gridded(temp_surf)
    u_grid = zero(temp_surf_grid)
    v_grid = zero(temp_surf_grid)
    vor_grid = zero(temp_surf_grid)

    P = Parameters(NF=Float32)
    G = GeoSpectral(P)

    SpeedyWeather.divergence_uvω_spectral(u_grid,v_grid,vor_grid,G)

end