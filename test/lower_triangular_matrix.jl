import JLArrays

@testset "LowerTriangularMatrix" begin
    @testset for NF in (Float32, Float64)
        mmax = 32
        @testset for lmax = (mmax, mmax+1)
            A = randn(Complex{NF}, lmax, mmax)

            SpeedyWeather.spectral_truncation!(A)

            L = SpeedyWeather.LowerTriangularMatrix(A)

            @test size(L) == size(A)

            for m in 1:mmax
                for l in 1:lmax
                    @test A[l, m] == L[l, m]
                end
            end

            @test_throws BoundsError L[lmax+1, mmax]
            @test_throws BoundsError L[lmax, mmax+1]
            
            @test_throws BoundsError L[1, 2] = 1
            @test_throws BoundsError L[lmax*mmax+1] = 1

            @test Matrix(L) == A
        end
    end
end

NF = Float32
mmax = 32 
idims = (5,)
lmax = mmax 
@testset "LowerTriangularArray: N-dim" begin 
    @testset for NF in (Float32, Float64)
        mmax = 32
        @testset for idims = ((), (5,), (5,5))
            @testset for lmax = (mmax, mmax+1)
                A = randn(Complex{NF}, lmax, mmax, idims...)

                # replaces spectraltrunction! here, we just set the elements zero manually
                ind = SpeedyWeather.LowerTriangularMatrices.lowertriangle_indices(A)
                ind = @. ~(ind)
                A[ind] .= zero(Complex{NF})

                L = SpeedyWeather.LowerTriangularArray(A)

                @test size(L) == size(A)

                for m in 1:mmax
                    for l in 1:lmax
                        @test A[l, m, [Colon() for i=1:length(idims)]...] == L[l, m, [Colon() for i=1:length(idims)]...]
                    end
                end

                @test_throws BoundsError L[lmax+1, mmax, [1 for i=1:length(idims)]...]
                @test_throws BoundsError L[lmax, mmax+1, [1 for i=1:length(idims)]...]
                
                @test_throws BoundsError L[1, 2, [1 for i=1:length(idims)]...] = 1
                @test_throws BoundsError L[lmax*mmax+1, [1 for i=1:length(idims)]...] = 1

                @test Array(L) == A
            end
        end
    end

end 

@testset "LowerTriangularMatrix: @inbounds" begin
    A = randn(LowerTriangularMatrix, 33, 32)
    
    @testset "getindex" begin
        @test_throws BoundsError A[34, 32]   # outside of i, j range
        @test_throws BoundsError A[561]     # outside of k range

        # # shouldn't throw an error if @inbounds
        # f(A, i) = @inbounds A[i]             # wrap into function
        # f(A, i, j) = @inbounds A[i, j]
        
        # f(A, 33, 32)                          # inside ranges
        # f(A, 560)

        # f(A, 34, 32)                          # outside ranges
        # f(A, 561)
    end

    @testset "setindex!" begin
        @test_throws BoundsError A[34, 32] = 1
        @test_throws BoundsError A[561] = 1

        # this can create a segfault so don't test regularly, but it worked ;)
        # g(A, i) = @inbounds A[i] = 1
        # g(A, i, j) = @inbounds A[i, j] = 1
        # g(A, 34, 32)
        # g(A, 561)
    end
end

@testset "LowerTriangularArray: @inbounds" begin
    A = randn(LowerTriangularArray{Float64}, 33, 32, 1, 1)
    
    @testset "getindex" begin
        @test_throws BoundsError A[34, 32, 1, 1]   # outside of i, j range
        @test_throws BoundsError A[561, 1, 1]     # outside of k range
        @test_throws BoundsError A[33, 32, 2, 1]     # outside of 3rd dim range
        @test_throws BoundsError A[33, 32, 1, 2]     # outside of 4th dim range
    end

    @testset "setindex!" begin
        @test_throws BoundsError A[34, 32, 1, 1] = 1
        @test_throws BoundsError A[561, 1, 1] = 1
        @test_throws BoundsError A[33, 32, 2, 1] = 1 
        @test_throws BoundsError A[33, 32, 1, 2] = 1 
    end
end

@testset "LowerTriangularArray: fill, copy, randn, convert" begin
    @testset for NF in (Float32, Float64)
        mmax = 32
        @testset for idims = ((), (5,), (5,5))
            @testset for lmax = (mmax, mmax+1)
                A = randn(Complex{NF}, lmax, mmax, idims...)

                # replaces spectraltrunction! here, we just set the elements zero manually
                ind = SpeedyWeather.LowerTriangularMatrices.lowertriangle_indices(A)
                ind = @. ~(ind)
                A[ind] .= zero(Complex{NF})

                L = SpeedyWeather.LowerTriangularArray(A)

                # fill
                fill!(L, 2)
                for lm in SpeedyWeather.eachharmonic(L)
                    @test all(L[lm, [Colon() for i=1:length(idims)]...] .== 2)
                end

                # copy
                L2 = copy(L)
                @test L2 == L

                L2[1, [1 for i=1:length(idims)]...] = 3
                @test L[1, [1 for i=1:length(idims)]...] == 2     # should be a deep copy
                @test L2[1, [1 for i=1:length(idims)]...] == 3

                # convert
                L = randn(LowerTriangularArray{NF}, lmax, mmax, idims...)
                L3 = convert(LowerTriangularArray{Float16, 2+length(idims), Array{Float16,1+length(idims)}}, L)
                for lm in SpeedyWeather.eachharmonic(L, L3)
                    @test Float16(L[lm, [1 for i=1:length(idims)]...]) == L3[lm, [1 for i=1:length(idims)]...] 
                end
            end
        end
    end
end

@testset "LowerTriangularMatrix: fill, copy, randn, convert" begin
    @testset for NF in (Float32, Float64)
        mmax = 32
        @testset for lmax = (mmax, mmax+1)
                A = randn(Complex{NF}, lmax, mmax)
                SpeedyWeather.spectral_truncation!(A)
                L = SpeedyWeather.LowerTriangularMatrix(A)

                # fill
                fill!(L, 2)
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
                L = randn(LowerTriangularMatrix{NF}, lmax, mmax)
                L3 = convert(LowerTriangularMatrix{Float16}, L)
                for lm in SpeedyWeather.eachharmonic(L, L3)
                    @test Float16(L[lm]) == L3[lm] 
                end
        end
    end
end

@testset "LowerTriangularMatrix: *, +, eachindex, similar" begin
    @testset for NF in (Float16, Float32, Float64)
        L = randn(LowerTriangularMatrix{NF}, 3, 3)

        @test (L+L) == 2L
        @test (L+L) == L*2
        @test (L-L) == zero(L)

        @test eachindex(L) == eachindex(L.data)
        @test eachindex(L, L) == eachindex(L.data, L.data)

        @test size(similar(L)) == size(L)
        @test eltype(L) == eltype(similar(L, eltype(L)))

        @test (5, 7) == size(similar(L, 5, 7))
        @test (5, 7) == size(similar(L, (5, 7)))
        @test similar(L) isa LowerTriangularMatrix
        @test similar(L, Float64) isa LowerTriangularMatrix{Float64}
    end
end

@testset "LowerTriangularArray: *, +, eachindex, similar" begin
    @testset for idims = ((), (5,), (5,5))
        @testset for NF in (Float16, Float32, Float64)
            L = randn(LowerTriangularArray{NF}, 3, 3, idims...)

            @test (L+L) == 2L
            @test (L+L) == L*2
            @test (L-L) == zero(L)

            @test eachindex(L) == eachindex(L.data)
            @test eachindex(L, L) == eachindex(L.data, L.data)

            @test size(similar(L)) == size(L)
            @test eltype(L) == eltype(similar(L, eltype(L)))

            @test (5, 7, idims...) == size(similar(L, 5, 7, idims...))
            @test (5, 7, idims...) == size(similar(L, (5, 7,  idims...)))
            @test similar(L) isa LowerTriangularArray
            @test similar(L, Float64) isa LowerTriangularArray{Float64}
        end
    end
end

@testset "LowerTriangularMatrix: copyto!" begin
    @testset for NF in (Float16, Float32, Float64)
        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = randn(LowerTriangularMatrix{NF}, 5, 5)
        L1c = copy(L1)

        copyto!(L2, L1)  # bigger into smaller
        copyto!(L1, L2)  # and back should be identical

        @test L1 == L1c

        # now smaller into bigger
        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = randn(LowerTriangularMatrix{NF}, 5, 5)
        L2c = copy(L2)

        copyto!(L1, L2)
        copyto!(L2, L1)

        @test L2 == L2c

        # with ranges
        L1 = zeros(LowerTriangularMatrix{NF}, 33, 32);
        L2 = randn(LowerTriangularMatrix{NF}, 65, 64);
        L2T = spectral_truncation(L2,(size(L1) .- 1)...)

        copyto!(L1, L2, 1:33, 1:32)     # size of smaller matrix
        @test L1 == L2T

        copyto!(L1, L2, 1:65, 1:64)     # size of bigger matrix
        @test L1 == L2T

        copyto!(L1, L2, 1:50, 1:50)     # in between
        @test L1 == L2T
    end
end
@testset "LowerTriangularArray: copyto!" begin
    @testset for idims = ((), (5,), (5,5))
        @testset for NF in (Float16, Float32, Float64)
            L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
            L2 = randn(LowerTriangularArray{NF}, 5, 5, idims...)
            L1c = copy(L1)

            copyto!(L2, L1)  # bigger into smaller
            copyto!(L1, L2)  # and back should be identical

            @test L1 == L1c

            # now smaller into bigger
            L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
            L2 = randn(LowerTriangularArray{NF}, 5, 5, idims...)
            L2c = copy(L2)

            copyto!(L1, L2)
            copyto!(L2, L1)

            @test L2 == L2c

            # with ranges
            L1 = zeros(LowerTriangularArray{NF}, 33, 32, idims...);
            L2 = randn(LowerTriangularArray{NF}, 65, 64, idims...);
            L2T = spectral_truncation(L2,(size(L1) .- 1)...)

            copyto!(L1, L2, 1:33, 1:32)     # size of smaller matrix
            @test L1 == L2T

            copyto!(L1, L2, 1:65, 1:64)     # size of bigger matrix
            @test L1 == L2T

            copyto!(L1, L2, 1:50, 1:50)     # in between
            @test L1 == L2T
        end
    end
end

@testset "LowerTriangularMatrix: broadcast" begin 
    @testset for NF in (Float16, Float32, Float64)
        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = copy(L1) 

        L2 .*= NF(5)
        @test L1 .* NF(5) ≈ L2 

        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = copy(L1) 

        L2 ./= NF(5)
        @test L1 ./ NF(5) ≈ L2 

        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = copy(L1)

        L2 .^= NF(2)
        @test L1 .^ NF(2) ≈ L2
    end
end 

@testset "LowerTriangularArray: broadcast" begin 
    @testset for idims = ((), (5,), (5,5))
        @testset for NF in (Float16, Float32, Float64)
            L1 = randn(LowerTriangularArray{NF, 2+length(idims), Array}, 10, 10, idims...)
            L2 = copy(L1) 

            L2 .*= NF(5)
            @test L1 .* NF(5) ≈ L2 

            L1 = randn(LowerTriangularArray{NF, 2+length(idims), Array}, 10, 10, idims...)
            L2 = copy(L1) 

            L2 ./= NF(5)
            @test L1 ./ NF(5) ≈ L2 

            L1 = randn(LowerTriangularArray{NF, 2+length(idims), Array}, 10, 10, idims...)
            L2 = copy(L1)

            L2 .^= NF(2)
            @test L1 .^ NF(2) ≈ L2
        end
    end
end 