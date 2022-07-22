@testset "Divergence of a non-divergent flow zero?" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater)

        fill!(p.vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.div,0)
        fill!(d.tendencies.vor_tend,0)

        # start with some vorticity only
        vor0 = randn(Complex{NF},size(p.vor)[1:2])
        SpeedyWeather.spectral_truncation!(vor0)
        p.vor[:,:,1,1] .= vor0       
        SpeedyWeather.gridded!(d,p,m)   # get corresponding non-divergent U_grid, V_grid

        U_grid = d.grid_variables.U_grid[:,:,1]
        V_grid = d.grid_variables.V_grid[:,:,1]

        # check we've actually created non-zero U=u*coslat,V=v*coslat
        @test all(U_grid .!= 0)
        @test all(V_grid .!= 0)

        G = m.geospectral.geometry
        S = m.geospectral.spectral_transform
        SpeedyWeather.scale_coslat⁻²!(U_grid,G)
        SpeedyWeather.scale_coslat⁻²!(V_grid,G)

        uω_coslat⁻¹ = d.intermediate_variables.uω_coslat⁻¹[:,:,1]
        vω_coslat⁻¹ = d.intermediate_variables.vω_coslat⁻¹[:,:,1]

        SpeedyWeather.spectral!(uω_coslat⁻¹,U_grid,S)
        SpeedyWeather.spectral!(vω_coslat⁻¹,V_grid,S)
    
        div = zero(vor0)
        SpeedyWeather.divergence!(div,uω_coslat⁻¹,vω_coslat⁻¹,S)

        for div_lm in div
            @test abs(div_lm) < sqrt(eps(NF))
        end
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
        @test all(d.grid_variables.U_grid .!= 0)
        @test all(d.grid_variables.V_grid .!= 0)

        # to evaluate ∇×(uv) use curl of vorticity fluxes (=∇×(uv(ζ+f))) with ζ=1,f=0
        fill!(d.grid_variables.vor_grid,1)
        fill!(m.geospectral.geometry.f_coriolis,0)

        # calculate uω,vω in spectral space
        SpeedyWeather.vorticity_flux_divergence!(d,m.geospectral)
        SpeedyWeather.vorticity_flux_curl!(d,m.geospectral)

        for div_lm in d.tendencies.div_tend
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Scale, unscale coslat" begin
    @testset for NF in (Float32,Float64)
        p,d,m = initialize_speedy(NF)
        G = m.geospectral.geometry

        A = randn(NF,96,48)
        B = copy(A)
        SpeedyWeather.scale_coslat⁻¹!(A,G)
        SpeedyWeather.scale_coslat!(A,G)

        @test all(isapprox.(A,B,rtol=10*eps(NF)))

        SpeedyWeather.scale_coslat²!(A,G)
        SpeedyWeather.scale_coslat⁻²!(A,G)

        @test all(isapprox.(A,B,rtol=10*eps(NF)))
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
                                    model=:barotropic)

        S = m.geospectral.spectral_transform
        A = randn(Complex{NF},size(p.vor)[1:2])
        B = zero(A)
        SpeedyWeather.gradient_latitude!(B,A,S,flipsign=true)
        @test SpeedyWeather.gradient_latitude(A,S,one_more_l=false) == -B
    end
end

@testset "Flipsign in divergence!, curl!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:barotropic)

        S = m.geospectral.spectral_transform
        lmax,mmax = size(p.vor)[1:2] .- 1
        A1 = randn(Complex{NF},lmax+2,mmax+1)
        A2 = randn(Complex{NF},lmax+2,mmax+1)
        B = zeros(Complex{NF},lmax+1,mmax+1)
        C = zeros(Complex{NF},lmax+1,mmax+1)

        SpeedyWeather.divergence!(B,A1,A2,S,flipsign=true)
        SpeedyWeather.divergence!(C,A1,A2,S,flipsign=false)
        @test C == -B

        SpeedyWeather.curl!(B,A1,A2,S,flipsign=true)
        SpeedyWeather.curl!(C,A1,A2,S,flipsign=false)
        @test C == -B
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

@testset "Add in divergence!, curl!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:barotropic)

        S = m.geospectral.spectral_transform
        lmax,mmax = size(p.vor)[1:2] .- 1
        A1 = randn(Complex{NF},lmax+2,mmax+1)
        A2 = randn(Complex{NF},lmax+2,mmax+1)
        B = zeros(Complex{NF},lmax+1,mmax+1)
        C = zeros(Complex{NF},lmax+1,mmax+1)

        SpeedyWeather.divergence!(B,A1,A2,S,add=true)
        SpeedyWeather.divergence!(B,A1,A2,S,add=true)
        SpeedyWeather.divergence!(C,A1,A2,S,add=false)
        @test 2C == B

        SpeedyWeather.divergence!(B,A1,A2,S,add=true)
        SpeedyWeather.divergence!(B,A1,A2,S,add=true,flipsign=true)
        @test all(2C .≈ B)

        fill!(B,0)
        SpeedyWeather.curl!(B,A1,A2,S,add=true)
        SpeedyWeather.curl!(B,A1,A2,S,add=true)
        SpeedyWeather.curl!(C,A1,A2,S,add=false)
        @test 2C == B

        SpeedyWeather.curl!(B,A1,A2,S,add=true)
        SpeedyWeather.curl!(B,A1,A2,S,add=true,flipsign=true)
        @test all(2C .≈ B)
    end
end

@testset "curl_div! same as curl! and divergence!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(NF)

        lmax,mmax,_,_ = size(p.vor)
        
        U = randn(Complex{NF},lmax+1,mmax)
        V = randn(Complex{NF},lmax+1,mmax)

        div0 = zeros(Complex{NF},lmax,mmax)
        div1 = zeros(Complex{NF},lmax,mmax)
        vor0 = zeros(Complex{NF},lmax,mmax)
        vor1 = zeros(Complex{NF},lmax,mmax)

        S = m.geospectral.spectral_transform
        SpeedyWeather.divergence!(div0,U,V,S)
        SpeedyWeather.curl!(vor0,U,V,S)

        SpeedyWeather.curl_div!(vor1,div1,U,V,S)

        for i in eachindex(vor0,vor1,div0,div1)
            @test vor0[i] == vor1[i]
            @test div0[i] == div1[i]
        end
    end
end

@testset "D,ζ -> u,v -> D,ζ" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater)

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

        vor1 = zero(vor0)
        div1 = zero(div0)

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        SpeedyWeather.gridded!(d,p,m)   

        # check we've actually created non-zero u,v
        @test all(d.grid_variables.U_grid .!= 0)
        @test all(d.grid_variables.V_grid .!= 0)

        U = d.grid_variables.U_grid[:,:,1]
        V = d.grid_variables.V_grid[:,:,1]

        # create coslat²*(div,vor) in spectral space
        div_grid = d.grid_variables.div_grid[:,:,1]
        vor_grid = d.grid_variables.vor_grid[:,:,1]

        # times coslat² in grid space
        G = m.geospectral.geometry
        SpeedyWeather.scale_coslat⁻²!(U,G)
        SpeedyWeather.scale_coslat⁻²!(V,G)

        # transform back
        S = m.geospectral.spectral_transform
        u_coslat⁻¹ = zero(d.intermediate_variables.u_coslat[:,:,1])
        v_coslat⁻¹ = zero(d.intermediate_variables.v_coslat[:,:,1])
        SpeedyWeather.spectral!(u_coslat⁻¹,U,S)
        SpeedyWeather.spectral!(v_coslat⁻¹,V,S)

        # curl and div in spectral space
        SpeedyWeather.curl!(vor1,u_coslat⁻¹,v_coslat⁻¹,S)
        SpeedyWeather.divergence!(div1,u_coslat⁻¹,v_coslat⁻¹,S)

        for i in eachindex(vor0,vor1,div0,div1)
            @test vor0[i] ≈ vor1[i] rtol=sqrt(eps(NF))
            @test div0[i] ≈ div1[i] rtol=sqrt(eps(NF))
        end
    end
end