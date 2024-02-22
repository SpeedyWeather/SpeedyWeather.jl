@testset "Divergence of a non-divergent flow zero?" begin
    @testset for NF in (Float32,Float64)

        spectral_grid = SpectralGrid(;NF,nlev=1)
        m = ShallowWaterModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        fill!(p.layers[1].timesteps[1].vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.layers[1].timesteps[1].div,0)
        fill!(d.layers[1].tendencies.vor_tend,0)

        # start with some vorticity only
        vor0 = randn(LowerTriangularMatrix{Complex{NF}},p.trunc+2,p.trunc+1)
        p.layers[1].timesteps[1].vor .= vor0

        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   # get corresponding non-divergent u_grid, v_grid

        u_grid = d.layers[1].grid_variables.u_grid
        v_grid = d.layers[1].grid_variables.v_grid

        # check we've actually created non-zero U=u*coslat,V=v*coslat
        @test all(u_grid .!= 0)
        @test all(v_grid .!= 0)

        G = m.geometry
        S = m.spectral_transform
        RingGrids.scale_coslat⁻¹!(u_grid)
        RingGrids.scale_coslat⁻¹!(v_grid)

        uω_coslat⁻¹ = d.layers[1].dynamics_variables.a
        vω_coslat⁻¹ = d.layers[1].dynamics_variables.b

        SpeedyWeather.spectral!(uω_coslat⁻¹,u_grid,S)
        SpeedyWeather.spectral!(vω_coslat⁻¹,v_grid,S)
    
        div = zero(vor0)
        SpeedyWeather.divergence!(div,uω_coslat⁻¹,vω_coslat⁻¹,S)

        for div_lm in div
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Curl of an irrotational flow zero?" begin
    @testset for NF in (Float32,Float64)

        spectral_grid = SpectralGrid(;NF,nlev=1)
        m = ShallowWaterModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        fill!(p.layers[1].timesteps[1].vor,0)                  # make sure vorticity and divergence are 0
        fill!(p.layers[1].timesteps[1].div,0)
        fill!(d.layers[1].tendencies.div_tend,0)

        # start with some vorticity only
        div0 = randn(LowerTriangularMatrix{Complex{NF}},p.trunc+2,p.trunc+1)
        p.layers[1].timesteps[1].div .= div0

        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   # get corresponding non-divergent u_grid, v_grid

        u_grid = d.layers[1].grid_variables.u_grid
        v_grid = d.layers[1].grid_variables.v_grid

        # check we've actually created non-zero U=u*coslat,V=v*coslat
        @test all(u_grid .!= 0)
        @test all(v_grid .!= 0)

        G = m.geometry
        S = m.spectral_transform
        C = m.constants

        # to evaluate ∇×(uv) use curl of vorticity fluxes (=∇×(uv(ζ+f))) with ζ=1,f=0
        fill!(d.layers[1].grid_variables.vor_grid,1)
        fill!(C.f_coriolis,0)

        # calculate uω,vω in spectral space
        SpeedyWeather.vorticity_flux_curldiv!(d.layers[1],C,G,S,div=true)

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

            SG = SpectralGrid(NF;Grid,nlev=1)
            G = Geometry(SG)

            A = Grid(randn(NF,SG.npoints))
            B = copy(A)
            RingGrids.scale_coslat⁻¹!(A)
            RingGrids.scale_coslat!(A)

            @test all(isapprox.(A,B,rtol=10*eps(NF)))

            RingGrids.scale_coslat²!(A)
            RingGrids.scale_coslat⁻²!(A)

            @test all(isapprox.(A,B,rtol=10*eps(NF)))
        end
    end
end

@testset "Flipsign in divergence!, curl!" begin
    @testset for NF in (Float32,Float64)

        SG = SpectralGrid(NF)
        S = SpectralTransform(SG)

        lmax,mmax = S.lmax,S.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)

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

        SG = SpectralGrid(NF)
        S = SpectralTransform(SG)

        lmax,mmax = S.lmax,S.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)

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

        spectral_grid = SpectralGrid(;NF,nlev=1)
        m = ShallowWaterModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        # make sure vorticity and divergence are 0
        fill!(p.layers[1].timesteps[1].vor,0)
        fill!(p.layers[1].timesteps[1].div,0)

        # make sure vorticity and divergence are 0
        fill!(d.layers[1].tendencies.vor_tend,0)                  
        fill!(d.layers[1].tendencies.div_tend,0)

        # create initial conditions
        (;lmax,mmax) = m.spectral_transform
        vor0 = rand(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        div0 = rand(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax+1)
        
        vor0[1,1] = 0                   # zero mean
        div0[1,1] = 0
        vor0[:,1] .= real(vor0[:,1])    # set imaginary component of m=0 to 0
        div0[:,1] .= real(div0[:,1])    # as the rotation of zonal modes is arbitrary
        vor0[end,:] .= 0                # set unusued last row (l=lmax+1) to zero
        div0[end,:] .= 0

        # copy into prognostic variables
        p.layers[1].timesteps[1].vor .= vor0
        p.layers[1].timesteps[1].div .= div0

        vor1 = zero(vor0)
        div1 = zero(div0)

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        lf = 1
        SpeedyWeather.gridded!(d,p,lf,m)   

        # check we've actually created non-zero u,v
        @test all(d.layers[1].grid_variables.u_grid .!= 0)
        @test all(d.layers[1].grid_variables.v_grid .!= 0)

        u = d.layers[1].grid_variables.u_grid
        v = d.layers[1].grid_variables.v_grid

        # times coslat² in grid space
        G = m.geometry
        RingGrids.scale_coslat⁻¹!(u)
        RingGrids.scale_coslat⁻¹!(v)

        # transform back
        S = m.spectral_transform
        u_coslat⁻¹ = zero(vor0)
        v_coslat⁻¹ = zero(vor0)
        SpeedyWeather.spectral!(u_coslat⁻¹,u,S)
        SpeedyWeather.spectral!(v_coslat⁻¹,v,S)

        # curl and div in spectral space
        SpeedyWeather.curl!(vor1,u_coslat⁻¹,v_coslat⁻¹,S)
        SpeedyWeather.divergence!(div1,u_coslat⁻¹,v_coslat⁻¹,S)

        for lm in SpeedyWeather.eachharmonic(vor0,vor1,div0,div1)
            # increased to 20 as 10 caused single fails every now and then
            @test vor0[lm] ≈ vor1[lm] rtol=20*sqrt(eps(NF))
            @test div0[lm] ≈ div1[lm] rtol=20*sqrt(eps(NF))
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

@testset "∇×∇=0 and ∇⋅∇=∇²" begin
    for NF in (Float32,Float64)

        trunc = 31
        spectral_grid = SpectralGrid(;NF,trunc,Grid=FullGaussianGrid,nlev=1)
        m = ShallowWaterModel(;spectral_grid)
        simulation = initialize!(m)
        p = simulation.prognostic_variables
        d = simulation.diagnostic_variables

        a = randn(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
        SpeedyWeather.spectral_truncation!(a)
        a[:,1] .= real.(a[:,1])

        dadx = zero(a)
        dady = zero(a)
        SpeedyWeather.∇!(dadx,dady,a,m.spectral_transform)

        dadx_grid = gridded(dadx,m.spectral_transform)
        dady_grid = gridded(dady,m.spectral_transform)
        
        RingGrids.scale_coslat⁻²!(dadx_grid)
        RingGrids.scale_coslat⁻²!(dady_grid)

        SpeedyWeather.spectral!(dadx,dadx_grid,m.spectral_transform)
        SpeedyWeather.spectral!(dady,dady_grid,m.spectral_transform)

        # CURL(GRAD(A)) = 0
        ∇x∇a = zero(a)
        SpeedyWeather.curl!(∇x∇a,dadx,dady,m.spectral_transform)

        for lm in SpeedyWeather.eachharmonic(∇x∇a)
            @test ∇x∇a[lm] ≈ 0 atol=5*sqrt(eps(NF))
        end
        
        # DIV(GRAD(A)) = LAPLACE(A)
        div_∇a = zero(a)
        SpeedyWeather.divergence!(div_∇a,dadx,dady,m.spectral_transform)
        ∇²a = zero(a)
        SpeedyWeather.∇²!(∇²a,a,m.spectral_transform)

        for lm in SpeedyWeather.eachharmonic(div_∇a,∇²a)
            @test div_∇a[lm] ≈ ∇²a[lm] atol=5*sqrt(eps(NF)) rtol=5*sqrt(eps(NF))
        end
    end
end