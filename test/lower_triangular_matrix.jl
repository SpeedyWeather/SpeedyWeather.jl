@testset "LowerTriangularMatrix" begin
    @testset for NF in (Float32,Float64)
        mmax = 32
        @testset for lmax = (mmax,mmax+1)
            A = randn(Complex{NF},lmax,mmax)
            SpeedyWeather.spectral_truncation!(A)

            L = SpeedyWeather.LowerTriangularMatrix(A)

            for m in 1:mmax
                for l in 1:lmax
                    @test A[l,m] == L[l,m]
                end
            end

            @test_throws BoundsError L[lmax+1,mmax]
            @test_throws BoundsError L[lmax,mmax+1]

            @test Matrix(L) == A
        end
    end
end