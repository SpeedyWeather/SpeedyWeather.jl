# @testset "Spectral gradients no error test" begin
#     for NF in (Float32,Float64)
#         prog_vars,diag_vars,model_setup = initialize_speedy(NF)

#         temp_surf = prog_vars.temp[:,:,1,end]
#         SpeedyWeather.gradient_longitude(temp_surf)
#         SpeedyWeather.gradient_latitude(temp_surf)

#         SpeedyWeather.∇²(temp_surf)
#         SpeedyWeather.∇⁻²(temp_surf)

#         temp_surf_grid = gridded(temp_surf)
        
#         u_surf_grid = diag_vars.grid_variables.u_grid[:,:,end]
#         v_surf_grid = diag_vars.grid_variables.v_grid[:,:,end]
#         vor_surf_grid = diag_vars.grid_variables.vor_grid[:,:,end]

#         SpeedyWeather.divergence_uvω_spectral(  u_surf_grid,
#                                                 v_surf_grid,
#                                                 vor_surf_grid,
#                                                 model_setup.geospectral)
#     end
# end

# @testset "gradient_latitude, with/without precalculated ϵlms" begin
#     for NF in (Float32,Float64)
#         prog_vars,diag_vars,model_setup = initialize_speedy(NF)

#         (;coslat_u,stream_function) = diag_vars.intermediate_variables
#         ϵlms = model_setup.geospectral.spectral.ϵlms

#         # some random entries
#         lmax,mmax = 5,5
#         stream_function[1:lmax,1:mmax,:] = convert(NF,1e-6)*randn(Complex{NF},lmax,mmax,
#                                                     model_setup.parameters.nlev)
#         SpeedyWeather.spectral_truncation!(stream_function,lmax,mmax)   # set upper triangle to 0
#         radius_earth = 2

#         # calculates ϵlms upfront
#         SpeedyWeather.gradient_latitude!(coslat_u,stream_function,radius_earth)
#         coslat_u2 = copy(coslat_u)      # store for comparison

#         # uses ϵlms from SpectralTransform
#         SpeedyWeather.gradient_latitude!(coslat_u,stream_function,ϵlms,radius_earth)

#         @test coslat_u == coslat_u2
#     end
# end

# @testset "Spectral vorticity->stream function->u,v->vorticity" begin
    # for NF in (Float32,Float64)
        NF = Float64
        prog_vars,diag_vars,model_setup = initialize_speedy(NF,trunc=85)
        (;spectral,geometry) = model_setup.geospectral
        (;radius_earth) = model_setup.parameters

        vor = view(prog_vars.vor,:,:,1,:)    # use only one leapfrog index

        # some large scale initial conditions
        lmax,mmax = 50,50
        vor[2:lmax,2:mmax,:] = randn(Complex{NF},lmax-1,mmax-1,
                                        model_setup.parameters.nlev)
        SpeedyWeather.spectral_truncation!(vor,lmax,mmax)   # set upper triangle to 0

        SpeedyWeather.gridded!(diag_vars,prog_vars,model_setup)

        (;u_grid, v_grid) = diag_vars.grid_variables
        u = zero(vor)
        v = zero(vor)

        SpeedyWeather.scale_coslat!(v_grid,geometry)

        SpeedyWeather.spectral!(u,u_grid,spectral)
        SpeedyWeather.spectral!(v,v_grid,spectral)

        (;coslat_u,coslat_v) = diag_vars.intermediate_variables
        SpeedyWeather.gradient_longitude!(coslat_v,v,radius_earth)
        SpeedyWeather.gradient_latitude!(coslat_u,u,spectral,radius_earth)

        vor2 = coslat_v - coslat_u
        vor_grid2 = gridded(vor2[:,:,1])
        SpeedyWeather.unscale_coslat!(vor_grid2,geometry)

        SpeedyWeather.spectral!(view(vor2,:,:,1),vor_grid2,spectral)

        fig,(ax1,ax2,ax3) = subplots(3,1)
        q1 = ax1.imshow(abs.(vor[:,:,1]))
        colorbar(q1,ax=ax1)
        q2 = ax2.imshow(abs.(vor2[:,:,1]))
        colorbar(q2,ax=ax2)
        q3 = ax3.imshow(abs.(vor[:,:,1])-abs.(vor2[:,:,1]))
        colorbar(q3,ax=ax3)
        ax1.set_title("Vorticity, random",loc="left")
        ax1.set_title("a",loc="right")
        ax2.set_title("Vorticity->stream function->u,v->vorticity")
        ax2.set_title("b",loc="right")
        ax3.set_title("Error: a-b")

    # end
# end

# @testset "Scale/unscale grid variables with coslat" begin
#     for NF in (Float32,Float64)
#         prog_vars,diag_vars,model_setup = initialize_speedy(NF)

#         (;u_grid) = diag_vars.grid_variables
#         u_grid = randn(NF,size(u_grid)...)
#         u_grid2 = copy(u_grid)

#         # first unscale then scale with coslat
#         SpeedyWeather.unscale_coslat!(u_grid2,model_setup.geospectral.geometry)
#         SpeedyWeather.scale_coslat!(u_grid2,model_setup.geospectral.geometry)

#         @test all(u_grid .≈ u_grid2)

#         # first scale then unscale with coslat
#         SpeedyWeather.scale_coslat!(u_grid2,model_setup.geospectral.geometry)
#         SpeedyWeather.unscale_coslat!(u_grid2,model_setup.geospectral.geometry)

#         @test all(u_grid .≈ u_grid2)
#     end
# end