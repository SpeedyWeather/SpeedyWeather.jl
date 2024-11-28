@testset "Divergence of a non-divergent flow zero?" begin
    @testset for NF in (Float32, Float64)

        spectral_grid = SpectralGrid(; NF, nlayers=1)
        model = ShallowWaterModel(spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables

        # make sure vorticity and divergence are 0
        progn.vor[1] .= 0
        progn.div[1] .= 0
        diagn.tendencies.vor_tend .= 0

        # start with some vorticity only
        progn.vor[1] .= randn(Complex{NF}, size(progn.vor[1])...)

        # get corresponding non-divergent u_grid, v_grid
        lf = 1
        transform!(diagn, progn, lf, model)

        (; u_grid, v_grid) = diagn.grid

        # check we've actually created non-zero u,v excl coslat scaling
        @test all(u_grid .!= 0)
        @test all(v_grid .!= 0)

        RingGrids.scale_coslat⁻¹!(u_grid)
        RingGrids.scale_coslat⁻¹!(v_grid)
        
        uω_coslat⁻¹ = diagn.dynamics.a
        vω_coslat⁻¹ = diagn.dynamics.b
        
        S = model.spectral_transform
        SpeedyWeather.transform!(uω_coslat⁻¹, u_grid, S)
        SpeedyWeather.transform!(vω_coslat⁻¹, v_grid, S)
    
        div = progn.div[1]
        SpeedyWeather.divergence!(div, uω_coslat⁻¹, vω_coslat⁻¹, S)

        for div_lm in div
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Curl of an irrotational flow zero?" begin
    @testset for NF in (Float32, Float64)

        spectral_grid = SpectralGrid(; NF, nlayers=1)
        model = ShallowWaterModel(spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables

        # make sure vorticity and divergence are 0
        progn.vor[1] .= 0
        progn.div[1] .= 0
        diagn.tendencies.div_tend .= 0

        # start with some divergence only
        progn.div[1] .= randn(Complex{NF}, size(progn.div[1])...)

        # get corresponding non-divergent u_grid, v_grid
        lf = 1
        transform!(diagn, progn, lf, model)

        (; u_grid, v_grid) = diagn.grid

        # check we've actually created non-zero u,v (excl coslat scaling)
        @test all(u_grid .!= 0)
        @test all(v_grid .!= 0)

        # to evaluate ∇×(uv) use curl of vorticity fluxes (=∇×(uv(ζ+f))) with ζ=1, f=0
        progn.div[1] .= 1
        model.coriolis.f .= 0

        # calculate uω, vω in spectral space
        (; coriolis, geometry, spectral_transform) = model
        SpeedyWeather.vorticity_flux_curldiv!(diagn, coriolis, geometry, spectral_transform, div=true)

        for div_lm in diagn.tendencies.div_tend
            @test abs(div_lm) < sqrt(eps(NF))
        end
    end
end

@testset "Scale, unscale coslat" begin
    @testset for NF in (Float32, Float64)
        for Grid in (   FullGaussianGrid,
                        FullClenshawGrid,
                        OctahedralGaussianGrid,
                        OctahedralClenshawGrid,
                        HEALPixGrid)

            SG = SpectralGrid(; NF, Grid, nlayers=1)

            A = Grid(randn(NF, SG.npoints))
            B = copy(A)
            RingGrids.scale_coslat⁻¹!(A)
            RingGrids.scale_coslat!(A)

            @test all(isapprox.(A, B, rtol=10*eps(NF)))

            RingGrids.scale_coslat²!(A)
            RingGrids.scale_coslat⁻²!(A)

            @test all(isapprox.(A, B, rtol=10*eps(NF)))
        end
    end
end

@testset "Flipsign in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)

        SG = SpectralGrid(; NF)
        S = SpectralTransform(SG)

        lmax, mmax = S.lmax, S.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)

        SpeedyWeather.divergence!(B, A1, A2, S, flipsign=true)
        SpeedyWeather.divergence!(C, A1, A2, S, flipsign=false)
        @test C == -B

        SpeedyWeather.curl!(B, A1, A2, S, flipsign=true)
        SpeedyWeather.curl!(C, A1, A2, S, flipsign=false)
        @test C == -B
    end
end

@testset "Add in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)

        SG = SpectralGrid(; NF)
        S = SpectralTransform(SG)

        lmax, mmax = S.lmax, S.mmax
        A1 = randn(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        A2 = randn(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        B = zeros(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)
        C = zeros(LowerTriangularMatrix{Complex{NF}}, lmax+1, mmax+1)

        SpeedyWeather.divergence!(B, A1, A2, S, add=true)
        SpeedyWeather.divergence!(B, A1, A2, S, add=true)
        SpeedyWeather.divergence!(C, A1, A2, S, add=false)
        @test 2C == B

        SpeedyWeather.divergence!(B, A1, A2, S, add=true)
        SpeedyWeather.divergence!(B, A1, A2, S, add=true, flipsign=true)
        @test all(2C .≈ B)

        fill!(B, 0)
        SpeedyWeather.curl!(B, A1, A2, S, add=true)
        SpeedyWeather.curl!(B, A1, A2, S, add=true)
        SpeedyWeather.curl!(C, A1, A2, S, add=false)
        @test 2C == B

        SpeedyWeather.curl!(B, A1, A2, S, add=true)
        SpeedyWeather.curl!(B, A1, A2, S, add=true, flipsign=true)
        @test all(2C .≈ B)
    end
end

@testset "D, ζ -> u, v -> D, ζ" begin
    @testset for NF in (Float32, Float64)

        spectral_grid = SpectralGrid(; NF, nlayers=2)
        model = PrimitiveDryModel(spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables
        vor = progn.vor[1]
        div = progn.div[1]

        # make sure vorticity and divergence are 0
        vor .= 0
        div .= 0

        # create initial conditions
        vor .= rand(Complex{NF}, size(vor)...)
        div .= rand(Complex{NF}, size(div)...)
        
        for k in eachmatrix(vor, div)
            vor[1, k] = 0     # zero mean on every layer
            div[1, k] = 0
        end

        # set imaginary component of m=0 to 0 as the rotation of zonal modes is arbitrary
        SpeedyTransforms.zero_imaginary_zonal_modes!(vor)
        SpeedyTransforms.zero_imaginary_zonal_modes!(div)

        spectral_truncation!(vor)      # set unusued last row (l=lmax+1) to zero
        spectral_truncation!(div)

        # get corresponding u_grid, v_grid (excl *coslat scaling)
        lf = 1
        SpeedyWeather.transform!(diagn, progn, lf, model)

        # check we've actually created non-zero u, v
        (; u_grid, v_grid) = diagn.grid
        @test all(u_grid .!= 0)
        @test all(v_grid .!= 0)

        # times coslat² in grid space
        RingGrids.scale_coslat⁻¹!(u_grid)
        RingGrids.scale_coslat⁻¹!(v_grid)

        # transform back
        u_coslat⁻¹ = transform(u_grid, model.spectral_transform)
        v_coslat⁻¹ = transform(v_grid, model.spectral_transform)

        # curl and div in spectral space
        vor2 = curl(u_coslat⁻¹, v_coslat⁻¹, model.spectral_transform)
        div2 = divergence(u_coslat⁻¹, v_coslat⁻¹, model.spectral_transform)

        for lm in eachindex(vor, div, vor2, div2)
            # increased to 30 as 10, 20 caused single fails every now and then
            @test vor[lm] ≈ vor2[lm] rtol=30*sqrt(eps(NF))
            @test div[lm] ≈ div2[lm] rtol=30*sqrt(eps(NF))
        end
    end
end

@testset "(Inverse) Laplace operator" begin

    for NF in (Float32, Float64)
        alms = LowerTriangularMatrix(randn(Complex{NF}, 33, 32))
        alms2 = copy(alms)
        alms3 = copy(alms)

        S = SpectralTransform(alms)

        # ∇⁻²! same as inverse=true
        SpeedyWeather.∇²!(alms2, alms, S, inverse=true);
        SpeedyWeather.∇⁻²!(alms3, alms, S);
        @test alms2 == alms3

        # test add=true
        fill!(alms2, 0)
        SpeedyWeather.∇²!(alms2, alms, S, add=true);
        SpeedyWeather.∇²!(alms3, alms, S);
        @test alms2 == alms3

        # also for inverse
        fill!(alms2, 0)
        SpeedyWeather.∇⁻²!(alms2, alms, S, add=true);
        SpeedyWeather.∇⁻²!(alms3, alms, S);
        @test alms2 == alms3

        # test flipsign
        SpeedyWeather.∇²!(alms2, alms, S, flipsign=true);
        SpeedyWeather.∇²!(alms3, alms, S);
        @test alms2 == -alms3

        # also for inverse
        SpeedyWeather.∇⁻²!(alms2, alms, S, flipsign=true);
        SpeedyWeather.∇⁻²!(alms3, alms, S);
        @test alms2 == -alms3

        # test ∇²(∇⁻²) = 1
        alms[1] = 0     # remove 0-mode which is set to zero
        SpeedyWeather.∇²!(alms2, alms, S);
        SpeedyWeather.∇⁻²!(alms3, alms2, S);
        @test alms ≈ alms3

        # and ∇⁻²(∇²) = 1
        SpeedyWeather.∇⁻²!(alms2, alms, S);
        SpeedyWeather.∇²!(alms3, alms2, S);
        @test alms ≈ alms3
    end
end

@testset "∇×∇=0 and ∇⋅∇=∇²" begin
    for nlayers in (1, 2)
        for NF in (Float32, Float64)

            trunc = 31
            spectral_grid = SpectralGrid(; NF, trunc, Grid=FullGaussianGrid, nlayers)
            model = PrimitiveDryModel(spectral_grid)
            simulation = initialize!(model)
            progn = simulation.prognostic_variables
            diagn = simulation.diagnostic_variables

            a = randn(LowerTriangularArray{Complex{NF}}, trunc+2, trunc+1, nlayers)
            spectral_truncation!(a)
            SpeedyTransforms.zero_imaginary_zonal_modes!(a)

            dadx, dady = ∇(a, model.spectral_transform)
            dadx_grid = transform(dadx, model.spectral_transform)
            dady_grid = transform(dady, model.spectral_transform)
            
            RingGrids.scale_coslat⁻²!(dadx_grid)
            RingGrids.scale_coslat⁻²!(dady_grid)

            transform!(dadx, dadx_grid, model.spectral_transform)
            transform!(dady, dady_grid, model.spectral_transform)

            # CURL(GRAD(A)) = 0
            ∇x∇a = curl(dadx, dady, model.spectral_transform)

            for lm in eachharmonic(∇x∇a)
                @test ∇x∇a[lm] ≈ 0 atol=5*sqrt(eps(NF))
            end
            
            # DIV(GRAD(A)) = LAPLACE(A)
            ∇dot∇a = divergence(dadx, dady, model.spectral_transform)
            ∇²a = ∇²(a, model.spectral_transform)

            for lm in eachindex(∇dot∇a, ∇²a)
                @test ∇dot∇a[lm] ≈ ∇²a[lm] atol=5*sqrt(eps(NF)) rtol=5*sqrt(eps(NF))
            end
        end
    end
end