@testset "Divergence of a non-divergent flow zero?" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=ShallowWater)

        fill!(p.layers[1].leapfrog[1].vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.layers[1].leapfrog[1].div,0)
        fill!(d.layers[1].tendencies.vor_tend,0)

        # start with some vorticity only
        vor0 = randn(LowerTriangularMatrix{Complex{NF}},p.lmax+2,p.mmax+1)
        p.layers[1].leapfrog[1].vor .= vor0

        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   # get corresponding non-divergent U_grid, V_grid

        U_grid = d.layers[1].grid_variables.U_grid
        V_grid = d.layers[1].grid_variables.V_grid

        # check we've actually created non-zero U=u*coslat,V=v*coslat
        @test all(U_grid .!= 0)
        @test all(V_grid .!= 0)

        G = m.geometry
        S = m.spectral_transform
        SpeedyWeather.scale_coslat⁻²!(U_grid,G)
        SpeedyWeather.scale_coslat⁻²!(V_grid,G)

        uω_coslat⁻¹ = d.layers[1].dynamics_variables.uω_coslat⁻¹
        vω_coslat⁻¹ = d.layers[1].dynamics_variables.vω_coslat⁻¹

        SpeedyWeather.spectral!(uω_coslat⁻¹,U_grid,S)
        SpeedyWeather.spectral!(vω_coslat⁻¹,V_grid,S)
    
        div = zero(vor0)
        SpeedyWeather.divergence!(div,uω_coslat⁻¹,vω_coslat⁻¹,S)

        for div_lm in div
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Curl of an irrotational flow zero?" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=ShallowWater)

        fill!(p.layers[1].leapfrog[1].vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.layers[1].leapfrog[1].div,0)
        fill!(d.layers[1].tendencies.div_tend,0)

        # start with some vorticity only
        div0 = randn(LowerTriangularMatrix{Complex{NF}},p.lmax+2,p.mmax+1)
        p.layers[1].leapfrog[1].div .= div0

        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   # get corresponding non-divergent U_grid, V_grid

        U_grid = d.layers[1].grid_variables.U_grid
        V_grid = d.layers[1].grid_variables.V_grid

        # check we've actually created non-zero U=u*coslat,V=v*coslat
        @test all(U_grid .!= 0)
        @test all(V_grid .!= 0)

        G = m.geometry
        S = m.spectral_transform

        # to evaluate ∇×(uv) use curl of vorticity fluxes (=∇×(uv(ζ+f))) with ζ=1,f=0
        fill!(d.layers[1].grid_variables.vor_grid,1)
        fill!(G.f_coriolis,0)

        # calculate uω,vω in spectral space
        SpeedyWeather.vorticity_flux_divergence!(d.layers[1],G,S)
        SpeedyWeather.vorticity_flux_curl!(d.layers[1],S)

        for div_lm in d.layers[1].tendencies.div_tend
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Scale, unscale coslat" begin
    @testset for NF in (Float32,Float64)
        for Grid in (   FullGaussianGrid,
                        FullClenshawGrid,
                        OctahedralGaussianGrid,
                        OctahedralClenshawGrid,
                        HEALPixGrid)

            p,d,m = initialize_speedy(NF;Grid)
            G = m.geometry

            A = Grid(randn(NF,G.npoints))
            B = copy(A)
            SpeedyWeather.scale_coslat⁻¹!(A,G)
            SpeedyWeather.scale_coslat!(A,G)

            @test all(isapprox.(A,B,rtol=10*eps(NF)))

            SpeedyWeather.scale_coslat²!(A,G)
            SpeedyWeather.scale_coslat⁻²!(A,G)

            @test all(isapprox.(A,B,rtol=10*eps(NF)))
        end
    end
end

@testset "Flipsign in divergence!, curl!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=BarotropicModel)

        S = m.spectral_transform
        lmax,mmax = p.lmax,p.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

        SpeedyWeather.divergence!(B,A1,A2,S,flipsign=true)
        SpeedyWeather.divergence!(C,A1,A2,S,flipsign=false)
        @test C == -B

        SpeedyWeather.curl!(B,A1,A2,S,flipsign=true)
        SpeedyWeather.curl!(C,A1,A2,S,flipsign=false)
        @test C == -B
    end
end

@testset "Add in divergence!, curl!" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=BarotropicModel)

        S = m.spectral_transform
        lmax,mmax = p.lmax,p.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)

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

@testset "D,ζ -> u,v -> D,ζ" begin
    @testset for NF in (Float32,Float64)

        p,d,m = initialize_speedy(  NF,
                                    model=ShallowWater)

        # make sure vorticity and divergence are 0
        fill!(p.layers[1].leapfrog[1].vor,0)
        fill!(p.layers[1].leapfrog[1].div,0)

        # make sure vorticity and divergence are 0
        fill!(d.layers[1].tendencies.vor_tend,0)                  
        fill!(d.layers[1].tendencies.div_tend,0)

        # create initial conditions
        lmax,mmax = p.lmax,p.mmax
        vor0 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        div0 = randn(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
        
        vor0[1,1] = 0                   # zero mean
        div0[1,1] = 0
        vor0[:,1] .= real(vor0[:,1])    # set imaginary component of m=0 to 0
        div0[:,1] .= real(div0[:,1])    # as the rotation of zonal modes is arbitrary
        vor0[end,:] .= 0                # set unusued last row (l=lmax+1) to zero
        div0[end,:] .= 0

        # copy into prognostic variables
        p.layers[1].leapfrog[1].vor .= vor0
        p.layers[1].leapfrog[1].div .= div0

        vor1 = zero(vor0)
        div1 = zero(div0)

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   

        # check we've actually created non-zero u,v
        @test all(d.layers[1].grid_variables.U_grid .!= 0)
        @test all(d.layers[1].grid_variables.V_grid .!= 0)

        U = d.layers[1].grid_variables.U_grid
        V = d.layers[1].grid_variables.V_grid

        # # create coslat²*(div,vor) in spectral space
        # div_grid = d.grid_variables.div_grid[:,:,1]
        # vor_grid = d.grid_variables.vor_grid[:,:,1]

        # times coslat² in grid space
        G = m.geometry
        SpeedyWeather.scale_coslat⁻²!(U,G)
        SpeedyWeather.scale_coslat⁻²!(V,G)

        # transform back
        S = m.spectral_transform
        u_coslat⁻¹ = zero(vor0)
        v_coslat⁻¹ = zero(vor0)
        SpeedyWeather.spectral!(u_coslat⁻¹,U,S)
        SpeedyWeather.spectral!(v_coslat⁻¹,V,S)

        # curl and div in spectral space
        SpeedyWeather.curl!(vor1,u_coslat⁻¹,v_coslat⁻¹,S)
        SpeedyWeather.divergence!(div1,u_coslat⁻¹,v_coslat⁻¹,S)

        for lm in SpeedyWeather.eachharmonic(vor0,vor1,div0,div1)
            @test vor0[lm] ≈ vor1[lm] rtol=10*sqrt(eps(NF))
            @test div0[lm] ≈ div1[lm] rtol=10*sqrt(eps(NF))
        end
    end
end

@testset "(Inverse) Laplace operator" begin

    for NF in (Float32,Float64)
        alms = LowerTriangularMatrix(randn(Complex{NF},33,32))
        alms2 = copy(alms)
        alms3 = copy(alms)

        S = SpectralTransform(alms,recompute_legendre=false)

        # ∇⁻²! same as inverse=true
        SpeedyWeather.∇²!(alms2,alms,S,inverse=true);
        SpeedyWeather.∇⁻²!(alms3,alms,S);
        @test alms2 == alms3

        # test add=true
        fill!(alms2,0)
        SpeedyWeather.∇²!(alms2,alms,S,add=true);
        SpeedyWeather.∇²!(alms3,alms,S);
        @test alms2 == alms3

        # also for inverse
        fill!(alms2,0)
        SpeedyWeather.∇⁻²!(alms2,alms,S,add=true);
        SpeedyWeather.∇⁻²!(alms3,alms,S);
        @test alms2 == alms3

        # test flipsign
        SpeedyWeather.∇²!(alms2,alms,S,flipsign=true);
        SpeedyWeather.∇²!(alms3,alms,S);
        @test alms2 == -alms3

        # also for inverse
        SpeedyWeather.∇⁻²!(alms2,alms,S,flipsign=true);
        SpeedyWeather.∇⁻²!(alms3,alms,S);
        @test alms2 == -alms3

        # test ∇²(∇⁻²) = 1
        alms[1] = 0     # remove 0-mode which is set to zero
        SpeedyWeather.∇²!(alms2,alms,S);
        SpeedyWeather.∇⁻²!(alms3,alms2,S);
        @test alms ≈ alms3

        # and ∇⁻²(∇²) = 1
        SpeedyWeather.∇⁻²!(alms2,alms,S);
        SpeedyWeather.∇²!(alms3,alms2,S);
        @test alms ≈ alms3
    end
end

@testset "∇ no errors" begin
    for NF in (Float32,Float64)
        NF = Float64
        p,d,m = initialize_speedy(NF)

        a = rand(LowerTriangularMatrix{Complex{NF}},33,32)
        dadx = zero(a)
        dady = zero(a)
        SpeedyWeather.∇!(dadx,dady,a,m.spectral_transform)
        
        # coslat unscaling missing in the following
        # ∇²a_1 = zero(a)
        # ∇²a_2 = zero(a)

        # SpeedyWeather.divergence!(∇²a_1,dadx,dady,m.spectral_transform)
        # SpeedyWeather.∇²!(∇²a_2,a,m.spectral_transform)

        # for lm in SpeedyWeather.eachharmonic(∇²a_1,∇²a_2)
        #     @test ∇²a_1[lm] ≈ ∇²a_2[lm]
        # end
    end
end