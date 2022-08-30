@testset "LowerTriangularMatrix" begin
    @testset for NF in (Float32,Float64)
        mmax = 32
        @testset for lmax = (mmax,mmax+1)
            A = randn(Complex{NF},lmax,mmax)
            SpeedyWeather.spectral_truncation!(A)

            L = SpeedyWeather.LowerTriangularMatrix(A)

            @test size(L) == size(A)

            for m in 1:mmax
                for l in 1:lmax
                    @test A[l,m] == L[l,m]
                end
            end

            @test_throws BoundsError L[lmax+1,mmax]
            @test_throws BoundsError L[lmax,mmax+1]
            
            @test_throws BoundsError L[1,2] = 1
            @test_throws BoundsError L[lmax*mmax+1] = 1

            @test Matrix(L) == A
        end
    end
end

@testset "LowerTriangularMatrix: fill, copy, randn, convert" begin
    @testset for NF in (Float32,Float64)
        mmax = 32
        @testset for lmax = (mmax,mmax+1)
            A = randn(Complex{NF},lmax,mmax)
            SpeedyWeather.spectral_truncation!(A)
            L = SpeedyWeather.LowerTriangularMatrix(A)

            # fill
            fill!(L,2)
            for lm in SpeedyWeather.eachharmonic(L)
                @test L[lm] == 2
            end

            # copy
            L2 = copy(L)
            @test L2 == L

            L2[1] = 3
            @test L[1] == 2     # should be a deep copy
            @test L2[1] == 3

            # convert
            L = randn(LowerTriangularMatrix{NF},lmax,mmax)
            L3 = convert(LowerTriangularMatrix{Float16},L)
            for lm in SpeedyWeather.eachharmonic(L,L3)
                @test Float16(L[lm]) == L3[lm] 
            end
        end
    end
end

@testset "LowerTriangularMatrix: *,+,eachindex,similar" begin
    @testset for NF in (Float16,Float32,Float64)
        L = LowerTriangularMatrix{NF}(undef,3,3)

        @test (L+L) == 2L
        @test (L+L) == L*2
        @test (L-L) == zero(L)

        @test eachindex(L) == eachindex(L.v)
        @test eachindex(L,L) == eachindex(L.v,L.v)

        @test size(similar(L)) == size(L)
        @test eltype(L) == eltype(similar(L,eltype(L)))

        @test (5,7) == size(similar(L,5,7))
        @test (5,7) == size(similar(L,(5,7)))
        @test similar(L) isa LowerTriangularMatrix
        @test similar(L,Float64) isa LowerTriangularMatrix{Float64}
    end
end
