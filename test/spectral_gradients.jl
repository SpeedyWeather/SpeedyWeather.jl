@testset "Divergence of a non-divergent flow zero?" begin

    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:barotropic,
                                    initial_conditions=:rest)

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.vor_tend,0)

        # start with some vorticity only
        p.vor[:,:,1,1] .= randn(Complex{NF},size(p.vor)[1:2])  
        SpeedyWeather.spectral_truncation!(p.vor[:,:,1,1])            
        SpeedyWeather.gridded!(d,p,m)   # get corresponding non-divergent u_grid, v_grid

        # check we've actually created non-zero u*coslat,v*coslat
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
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0)

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.pres_tend,0) # and tendency where the result goes
        fill!(m.boundaries.orography,0) # no mountains

        # start with some vorticity only
        p.vor[:,:,1,1] .= randn(Complex{NF},size(p.vor)[1:2])  
        SpeedyWeather.spectral_truncation!(p.vor[:,:,1,1])      
        SpeedyWeather.gridded!(d,p,m)   # get corresponding non-divergent u_grid, v_grid

        # check we've actually created non-zero u,v (which include coslat scaling)
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
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0)

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.div_tend,0)  # and tendency where the result goes

        # start with some divergence only
        p.div[:,:,1,1] .= randn(Complex{NF},size(p.div)[1:2])  
        SpeedyWeather.spectral_truncation!(p.div[:,:,1,1])

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        SpeedyWeather.gridded!(d,p,m)   

        # check we've actually created non-zero (u,v)*coslat
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

@testset "Flipsign in gradient_longitude!" begin
    @testset for NF in (Float32,Float64)
        A = randn(Complex{NF},32,32)
        B = zero(A)
        SpeedyWeather.gradient_longitude!(B,A,flipsign=true)
        @test SpeedyWeather.gradient_longitude(A,one_more_l=false) == -B
    end
end

@testset "Flipsign in gradient_latitude!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater)

        S = m.geospectral.spectral_transform
        A = randn(Complex{NF},size(p.vor)[1:2])
        B = zero(A)
        SpeedyWeather.gradient_latitude!(B,A,S,flipsign=true)
        @test SpeedyWeather.gradient_latitude(A,S,one_more_l=false) == -B
    end
end

@testset "Add in gradient_longitude!" begin
    @testset for NF in (Float32,Float64)
        A = randn(Complex{NF},32,32)
        SpeedyWeather.spectral_truncation!(A)

        B = zero(A)
        B2 = zero(A)
        SpeedyWeather.gradient_longitude!(B,A,add=false)
        SpeedyWeather.gradient_longitude!(B2,A,add=true)
        @test B == B2       # starting from 0 add doesn't make a difference?

        SpeedyWeather.gradient_longitude!(B2,A,add=true)
        @test 2B == B2      # adding the gradient twice is =x2?

        SpeedyWeather.gradient_longitude!(B2,A,add=true,flipsign=true)
        @test B == B2       # subtracting the gradient again returns to 1x gradient?
    end
end

@testset "Add in gradient_latitude!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(NF)

        S = m.geospectral.spectral_transform

        A = randn(Complex{NF},size(p.vor)[1:2])
        SpeedyWeather.spectral_truncation!(A)

        B = zero(A)
        B2 = zero(A)
        SpeedyWeather.gradient_latitude!(B,A,S,add=false)
        SpeedyWeather.gradient_latitude!(B2,A,S,add=true)
        @test B == B2       # starting from 0 add doesn't make a difference?

        SpeedyWeather.gradient_latitude!(B2,A,S,add=true)
        @test 2B == B2      # adding the gradient twice is =x2?

        SpeedyWeather.gradient_latitude!(B2,A,S,add=true,flipsign=true)
        @test B == B2       # subtracting the gradient again returns to 1x gradient?
    end
end

@testset "D,ζ -> u,v -> D,ζ" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0)

        # make sure vorticity and divergence are 0
        fill!(p.vor,0)
        fill!(p.div,0)

        # make sure vorticity and divergence are 0
        fill!(d.tendencies.vor_tend,0)                  
        fill!(d.tendencies.div_tend,0)

        # create initial conditions
        vor0 = zero(p.vor[:,:,1,1])
        div0 = zero(p.div[:,:,1,1])
        vor0[:,:] .= randn(Complex{NF},size(vor0)...)
        div0[:,:] .= randn(Complex{NF},size(div0)...)
        
        vor0[1,1] = 0                   # zero mean
        div0[1,1] = 0
        vor0[:,1] .= real(vor0[:,1])    # set imaginary component of m=0 to 0
        div0[:,1] .= real(div0[:,1])    # as the rotation of zonal modes is arbitrary

        # remove non-zero entries in upper triangle
        SpeedyWeather.spectral_truncation!(vor0)
        SpeedyWeather.spectral_truncation!(div0)

        p.vor[:,:,1,1] .= vor0
        p.div[:,:,1,1] .= div0

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        SpeedyWeather.gridded!(d,p,m)   

        # check we've actually created non-zero u,v
        @test all(d.grid_variables.u_grid .!= 0)
        @test all(d.grid_variables.v_grid .!= 0)

        coslat_u = d.intermediate_variables.coslat_u[:,:,1]
        coslat_v = d.intermediate_variables.coslat_v[:,:,1]

        # check that highest degree is non-zero
        @test all(coslat_u[end,:] .!= 0)
        @test all(coslat_v[end,:] .!= 0)

        # create coslat²*(div,vor) in spectral space
        div_grid = d.grid_variables.div_grid[:,:,1]
        vor_grid = d.grid_variables.vor_grid[:,:,1]

        # times coslat² in grid space
        G = m.geospectral.geometry
        SpeedyWeather.scale_coslat!(div_grid,G)
        SpeedyWeather.scale_coslat!(div_grid,G)
        SpeedyWeather.scale_coslat!(vor_grid,G)
        SpeedyWeather.scale_coslat!(vor_grid,G)

        # transform back
        S = m.geospectral.spectral_transform
        coslat²_div = spectral(div_grid,S,one_more_l=true)
        coslat²_vor = spectral(vor_grid,S,one_more_l=true)

        # zonal derivative in spectral space
        dUdlon = SpeedyWeather.gradient_longitude(coslat_u)
        dVdlon = SpeedyWeather.gradient_longitude(coslat_v)

        # meridional derivative of U,V in spectral space
        coslat_dVdθ = SpeedyWeather.gradient_latitude(coslat_v,S)
        coslat_dUdθ = SpeedyWeather.gradient_latitude(coslat_u,S)

        # test that
        # 1) coslat_dV/dlat = coslat²*D - dU/dlon
        # 2) coslat_dU/dlat = dV/dlon - coslat²*ζ
        for i in eachindex( coslat²_div,coslat²_vor,
                            coslat_dVdθ,coslat_dUdθ,
                            dUdlon,dVdlon)
            @test coslat_dVdθ[i] ≈ (coslat²_div[i] - dUdlon[i])
            @test coslat_dUdθ[i] ≈ (dVdlon[i] - coslat²_vor[i]) 
        end

        # # PLOTTING
        # fig,axs = subplots(3,2,figsize=(8,9))
        # levs = (0,3)
        # axs[1,1].imshow(abs.(coslat_dVdθ),vmin=levs[1],vmax=levs[2])
        # axs[1,2].imshow(abs.(coslat_dUdθ),vmin=levs[1],vmax=levs[2])
        # axs[2,1].imshow(abs.(coslat²_div),vmin=levs[1],vmax=levs[2])
        # q1 = axs[2,2].imshow(abs.(coslat²_vor),vmin=levs[1],vmax=levs[2])

        # levs2 = (0,0.1)
        # axs[3,1].imshow(abs.(coslat_dVdθ-coslat²_div),vmin=levs2[1],vmax=levs2[2])
        # q2 = axs[3,2].imshow(abs.(coslat_dUdθ-coslat²_vor),vmin=levs2[1],vmax=levs2[2])
        # colorbar(ax=axs[1:2,:],q1)
        # colorbar(ax=axs[3,:],q2)
    end
end