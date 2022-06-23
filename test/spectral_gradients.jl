@testset "Divergence of a non-divergent flow zero?" begin

    for NF in (Float32,Float64)

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

    for NF in (Float32,Float64)

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

    for NF in (Float32,Float64)

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

@testset "Vor,Div -> u,v -> vor,div" begin
    @testset for NF in (Float32,)#Float64)

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

        # G = m.geospectral.geometry
        # # SpeedyWeather.unscale_coslat!(d.grid_variables.u_grid,G)
        # # SpeedyWeather.unscale_coslat!(d.grid_variables.u_grid,G)
        # # SpeedyWeather.unscale_coslat!(d.grid_variables.v_grid,G)
        # # SpeedyWeather.unscale_coslat!(d.grid_variables.v_grid,G)

        S = m.geospectral.spectral_transform
        SpeedyWeather.spectral!(d.intermediate_variables.coslat_u,d.grid_variables.u_grid,S)
        SpeedyWeather.spectral!(d.intermediate_variables.coslat_v,d.grid_variables.v_grid,S)

        # TEST VORTICITY
        SpeedyWeather.gradient_longitude!(d.tendencies.vor_tend,d.intermediate_variables.coslat_v)
        SpeedyWeather.gradient_latitude!( d.tendencies.vor_tend,d.intermediate_variables.coslat_u,S,flipsign=true,add=true)

        lmax,mmax = size(vor0) .- 1
        for m in 1:mmax+1
            for l in m:lmax+1
                @test d.tendencies.vor_tend[l,m] ≈ vor0[l,m]
            end
        end

        # TEST DIVERGENCE
        SpeedyWeather.gradient_longitude!(d.tendencies.div_tend,d.intermediate_variables.coslat_u)
        SpeedyWeather.gradient_latitude!( d.tendencies.div_tend,d.intermediate_variables.coslat_v,S,add=true)

        lmax,mmax = size(div0) .- 1
        for m in 1:mmax+1
            for l in m:lmax+1
                @test d.tendencies.div_tend[l,m] ≈ div0[l,m]
            end
        end
    end
end