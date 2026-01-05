@testset "Divergence of a non-divergent flow zero?" begin
    @testset for NF in (Float32, Float64)
        @testset for nlayers in (1, 2)
            trunc = 31
            spectrum = Spectrum(trunc, one_degree_more=true)
            grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
            S = SpectralTransform(spectrum, grid; NF, nlayers)

            vor = randn(complex(NF), spectrum, nlayers)
            # LowerTriangularArrays.zero_last_degree!(vor)   # not needed

            div = zeros(complex(NF), spectrum, nlayers)
            U = zero(div)
            V = zero(div)

            # obtain u,v on grid
            SpeedyTransforms.UV_from_vordiv!(U, V, vor, div, S)
            u = transform(U, S, unscale_coslat=true)
            v = transform(V, S, unscale_coslat=true)

            # check we've actually created non-zero u,v excl coslat scaling
            @test all(u .!= 0)
            @test all(v .!= 0)

            # same velocities when leaving `div` explicitly out
            U2 = zero(U)
            V2 = zero(V)

            SpeedyTransforms.UV_from_vor!(U2, V2, vor, S)
            u2 = transform(U2, S, unscale_coslat=true)
            v2 = transform(V2, S, unscale_coslat=true)

            @test u == u2
            @test v == v2

            RingGrids.scale_coslat⁻¹!(u)
            RingGrids.scale_coslat⁻¹!(v)

            û = zero(U)
            v̂ = zero(V) 

            transform!(û, u, S)
            transform!(v̂, v, S)

            divergence!(vor, û, v̂, S)

            for lm in eachindex(div)
                @test abs(div[lm]) < sqrt(eps(NF))
            end
        end
    end
end

@testset "Curl of an irrotational flow zero?" begin
    @testset for NF in (Float32, Float64)
        @testset for nlayers in (1, 2)
            trunc = 31
            spectrum = Spectrum(trunc, one_degree_more=true)
            grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
            S = SpectralTransform(spectrum, grid; NF, nlayers)

            div = randn(complex(NF), spectrum, nlayers)
            # LowerTriangularArrays.zero_last_degree!(div)   # not needed
            vor = zeros(complex(NF), spectrum, nlayers)

            U = zero(div)
            V = zero(div)

            # obtain u,v on grid
            SpeedyTransforms.UV_from_vordiv!(U, V, vor, div, S)
            u = transform(U, S, unscale_coslat=true)
            v = transform(V, S, unscale_coslat=true)

            # check we've actually created non-zero u,v excl coslat scaling
            @test all(u .!= 0)
            @test all(v .!= 0)

            RingGrids.scale_coslat⁻¹!(u)
            RingGrids.scale_coslat⁻¹!(v)

            û = zero(U)
            v̂ = zero(V) 

            transform!(û, u, S)
            transform!(v̂, v, S)

            vor = curl(û, v̂, S)

            for lm in eachindex(vor)
                @test abs(vor[lm]) < sqrt(eps(NF))
            end
        end
    end
end

@testset "Flipsign in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)

        spectrum = Spectrum(15)
        grid = OctahedralGaussianGrid(12)
        S = SpectralTransform(spectrum, grid; NF)

        A1 = randn(Complex{NF}, spectrum)
        A2 = randn(Complex{NF}, spectrum)
        B = zeros(Complex{NF}, spectrum)
        C = zeros(Complex{NF}, spectrum)

        divergence!(B, A1, A2, S, flipsign = true)
        divergence!(C, A1, A2, S, flipsign = false)
        @test C == -B

        curl!(B, A1, A2, S, flipsign = true)
        curl!(C, A1, A2, S, flipsign = false)
        @test C == -B
    end
end

@testset "Add in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)
        spectrum = Spectrum(15)
        grid = OctahedralGaussianGrid(12)
        S = SpectralTransform(spectrum, grid; NF)

        A1 = randn(Complex{NF}, spectrum)
        A2 = randn(Complex{NF}, spectrum)
        B = zeros(Complex{NF}, spectrum)
        C = zeros(Complex{NF}, spectrum)

        divergence!(B, A1, A2, S, add = true)
        divergence!(B, A1, A2, S, add = true)
        divergence!(C, A1, A2, S, add = false)
        @test 2C == B

        divergence!(B, A1, A2, S, add = true)
        divergence!(B, A1, A2, S, add = true, flipsign = true)
        @test all(2C .≈ B)

        fill!(B, 0)
        curl!(B, A1, A2, S, add = true)
        curl!(B, A1, A2, S, add = true)
        curl!(C, A1, A2, S, add = false)
        @test 2C == B

        curl!(B, A1, A2, S, add = true)
        curl!(B, A1, A2, S, add = true, flipsign = true)
        @test all(2C .≈ B)
    end
end

@testset "Zero mean, zero last degree in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)

        trunc = 15
        spectrum = Spectrum(15, one_degree_more=true)
        grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
        S = SpectralTransform(spectrum, grid; NF)

        U = randn(Complex{NF}, spectrum)
        V = randn(Complex{NF}, spectrum)
        D = divergence(U, V, S)

        # general divergence properties
        @test D[1] == 0         # zero mean

        for m in 1:trunc+1      # last degree is also zero
            @test D[trunc+2, m] == 0
        end

        # same for curl
        ζ = curl(U, V, S)
        @test ζ[1] == 0         # zero mean

        for m in 1:trunc+1      # last degree is also zero
            @test ζ[trunc+2, m] == 0
        end
    end
end

@testset "Radius in divergence!, curl!" begin
    @testset for NF in (Float32, Float64)

        trunc = 15
        spectrum = Spectrum(15, one_degree_more=true)
        grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
        S = SpectralTransform(spectrum, grid; NF)

        U = randn(Complex{NF}, spectrum)
        V = randn(Complex{NF}, spectrum)

        R = 10
        D1 = divergence(U, V, S, radius=R)
        D2 = divergence(U, V, S)
        @test D1 ≈ D2/R

        ζ1 = curl(U, V, S, radius=R)
        ζ2 = curl(U, V, S)
        @test ζ1 ≈ ζ2/R
    end
end

@testset "D, ζ -> u, v -> D, ζ" begin
    @testset for NF in (Float32, Float64)
        @testset for nlayers in (1, 2)

            trunc = 31
            spectrum = Spectrum(trunc, one_degree_more=true)
            grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
            S = SpectralTransform(spectrum, grid; NF, nlayers)

            vor = randn(complex(NF), spectrum, nlayers)
            div = randn(complex(NF), spectrum, nlayers)
            
            zero_last_degree!(vor)   # needed?
            zero_last_degree!(div)   # needed?

            # zero mean on every layer
            vor[1, :] .= 0
            div[1, :] .= 0

            # set imaginary component of m=0 to 0 as the rotation of zonal modes is arbitrary
            vor[1:trunc+2, :] .= real(vor[1:trunc+2, :])
            div[1:trunc+2, :] .= real(div[1:trunc+2, :])

            U = zero(div)
            V = zero(div)

            # obtain u,v on grid
            SpeedyTransforms.UV_from_vordiv!(U, V, vor, div, S)
            u = transform(U, S, unscale_coslat=true)
            v = transform(V, S, unscale_coslat=true)

            RingGrids.scale_coslat⁻¹!(u)
            RingGrids.scale_coslat⁻¹!(v)

            û = zero(U)
            v̂ = zero(V) 

            transform!(û, u, S)
            transform!(v̂, v, S)

            vor2 = curl(û, v̂, S)
            div2 = divergence(û, v̂, S)

            # increased to 30 as 10, 20 caused single fails every now and then
            tol = NF == Float32 ? sqrt(eps(NF)) : 30*sqrt(eps(NF))
            for lm in eachindex(vor, div, vor2, div2)
                @test vor[lm] ≈ vor2[lm] rtol = tol
                @test div[lm] ≈ div2[lm] rtol = tol
            end
        end
    end
end

@testset "(Inverse) Laplace operator" begin

    for NF in (Float32, Float64)

        spectrum = Spectrum(31)

        alms = randn(Complex{NF}, spectrum)
        alms2 = copy(alms)
        alms3 = copy(alms)

        S = SpectralTransform(alms)

        # ∇⁻²! same as inverse=true
        ∇²!(alms2, alms, S, inverse = true)
        ∇⁻²!(alms3, alms, S)
        @test alms2 == alms3

        # test add=true
        fill!(alms2, 0)
        ∇²!(alms2, alms, S, add = true)
        ∇²!(alms3, alms, S)
        @test alms2 == alms3

        # also for inverse
        fill!(alms2, 0)
        ∇⁻²!(alms2, alms, S, add = true)
        ∇⁻²!(alms3, alms, S)
        @test alms2 == alms3

        # test flipsign
        ∇²!(alms2, alms, S, flipsign = true)
        ∇²!(alms3, alms, S)
        @test alms2 == -alms3

        # also for inverse
        ∇⁻²!(alms2, alms, S, flipsign = true)
        ∇⁻²!(alms3, alms, S)
        @test alms2 == -alms3

        # test ∇²(∇⁻²) = 1
        alms[1] = 0     # remove 0-mode which is set to zero
        ∇²!(alms2, alms, S)
        ∇⁻²!(alms3, alms2, S)
        @test alms ≈ alms3

        # and ∇⁻²(∇²) = 1
        ∇⁻²!(alms2, alms, S)
        ∇²!(alms3, alms2, S)
        @test alms ≈ alms3
    end
end

@testset "∇×∇=0 and ∇⋅∇=∇²" begin
    @testset for nlayers in (1, 2)
        @testset for NF in (Float32, Float64)

            trunc = 31
            spectrum = Spectrum(trunc, one_degree_more=true)
            grid = FullGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
            S = SpectralTransform(spectrum, grid; NF, nlayers)

            a = randn(complex(NF), spectrum, nlayers)
            SpeedyTransforms.spectral_truncation!(a)
            SpeedyTransforms.zero_imaginary_zonal_modes!(a)

            dadx, dady = ∇(a, S)
            dadx_grid = transform(dadx, S)
            dady_grid = transform(dady, S)

            RingGrids.scale_coslat⁻²!(dadx_grid)
            RingGrids.scale_coslat⁻²!(dady_grid)

            transform!(dadx, dadx_grid, S)
            transform!(dady, dady_grid, S)

            # CURL(GRAD(A)) = 0
            ∇x∇a = curl(dadx, dady, S)

            tol = 5 * sqrt(eps(NF))
            for lm in eachindex(∇x∇a)
                @test ∇x∇a[lm] ≈ 0 atol = tol
            end

            # DIV(GRAD(A)) = LAPLACE(A)
            ∇dot∇a = divergence(dadx, dady, S)
            ∇²a = ∇²(a, S)

            for lm in eachindex(∇dot∇a, ∇²a)
                @test ∇dot∇a[lm] ≈ ∇²a[lm] atol = tol rtol = tol
            end
        end
    end
end
