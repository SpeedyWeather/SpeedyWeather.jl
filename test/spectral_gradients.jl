@testset "Divergence of a non-divergent flow zero?" begin

    for NF in (Float32,Float64)

        p,d,m = initialize_speedy(NF,initial_conditions=:rest);

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.vor_tend,0)

        p.vor[5,4,1,1] = 1              # start with some vorticity only
        SpeedyWeather.gridded!(d,p,m)   # get corresponding non-divergent u_grid, v_grid

        # check we've actually created non-zero u,v
        @test all(d.grid_variables.u_grid .!= 0)
        @test all(d.grid_variables.v_grid .!= 0)

        # to evaluate ∇⋅(uv) use vorticity adv (=∇⋅(uv(ζ+f))) with ζ=1,f=0
        fill!(d.grid_variables.vor_grid,1)
        fill!(m.geospectral.geometry.f_coriolis,0)
        SpeedyWeather.vorticity_advection!(d,m.geospectral)

        @test all(abs.(d.tendencies.vor_tend) .< sqrt(eps(NF)))
    end
end

@testset "Continuity: Div of a non-divergent flow zero?" begin

    for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0);

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.pres_tend,0) # and tendency where the result goes
        fill!(m.boundaries.orography,0) # no mountains

        p.vor[5,4,1,1] = 1              # start with some vorticity only
        SpeedyWeather.gridded!(d,p,m)   # get corresponding non-divergent u_grid, v_grid

        # check we've actually created non-zero u,v
        @test all(d.grid_variables.u_grid .!= 0)
        @test all(d.grid_variables.v_grid .!= 0)
        
        # to evaluate ∇⋅(uv) use vorticity adv (=∇⋅(uv(η+H₀))) with η=1,H₀=0
        H₀ = 0
        fill!(d.grid_variables.pres_grid,1)
        SpeedyWeather.volume_fluxes!(d,m.geospectral,m.boundaries,H₀)

        @test all(abs.(d.tendencies.pres_tend) .< sqrt(eps(NF)))
    end
end

@testset "Curl of a irrotational flow zero?" begin

    for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0);

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.div_tend,0)  # and tendency where the result goes

        p.div[5,4,1,1] = 1              # start with some divergence only
        SpeedyWeather.gridded!(d,p,m)   # get corresponding irrotational u_grid, v_grid

        # check we've actually created non-zero u,v
        @test all(d.grid_variables.u_grid .!= 0)
        @test all(d.grid_variables.v_grid .!= 0)

        # to evaluate ∇×(uv) use curl of vorticity fluxes (=∇×(uv(ζ+f))) with ζ=1,f=0
        fill!(d.grid_variables.vor_grid,1)
        fill!(m.geospectral.geometry.f_coriolis,0)

        # calculate uω,vω in spectral space
        SpeedyWeather.vorticity_advection!(d,m.geospectral)
        SpeedyWeather.curl_vorticity_fluxes!(d,m.geospectral)

        @test all(abs.(d.tendencies.div_tend) .< sqrt(eps(NF)))
    end
end