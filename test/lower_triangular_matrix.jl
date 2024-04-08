using JLArrays
using Adapt

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

            @test all(L[:,1] .== A[:,1])
            @test all(L[:,2] .== A[:,2])
        end
    end
end

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

                if length(idims) == 1 
                    L1 = L[:,1] 
                    @test typeof(L1) <: LowerTriangularArray
                    @test all(L1.data .== L.data[:,1])
                    @test all([L1[i] == L[i,1] for i=1:size(L.data,1)])
                end
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
                L2 = deepcopy(L)
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
        L1c = deepcopy(L1)

        copyto!(L2, L1)  # bigger into smaller
        copyto!(L1, L2)  # and back should be identical

        @test L1 == L1c

        # now smaller into bigger
        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = randn(LowerTriangularMatrix{NF}, 5, 5)
        L2c = deepcopy(L2)

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

            # copyto! same size 
            L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
            L2 = similar(L1)
            copyto!(L2, L1)

            @test all(L2 .== L1)
            
            # copyto! Array 
            M = zeros(NF, 10, 10, idims...)
            copyto!(M, L1)

            ind = SpeedyWeather.LowerTriangularMatrices.lowertriangle_indices(M)
            not_ind = @. ~(ind)

            @test all(M[not_ind] .== zero(NF))
            @test all(LowerTriangularArray(M) .== L1)

            L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
            L2 = randn(LowerTriangularArray{NF}, 5, 5, idims...)
            L1c = deepcopy(L1)

            copyto!(L2, L1)  # bigger into smaller
            copyto!(L1, L2)  # and back should be identical

            @test L1 == L1c

            # now smaller into bigger
            L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
            L2 = randn(LowerTriangularArray{NF}, 5, 5, idims...)
            L2c = deepcopy(L2)

            copyto!(L1, L2)
            copyto!(L2, L1)

            @test L2 == L2c

            # with ranges
            L1 = zeros(LowerTriangularArray{NF}, 33, 32, idims...);
            L2 = randn(LowerTriangularArray{NF}, 65, 64, idims...);
            L2T = spectral_truncation(L2,(size(L1)[1:2] .- 1)...)

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
        L2 = deepcopy(L1) 

        L2 .*= NF(5)
        @test L1 .* NF(5) ≈ L2 

        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = deepcopy(L1) 

        L2 ./= NF(5)
        @test L1 ./ NF(5) ≈ L2 

        L1 = randn(LowerTriangularMatrix{NF}, 10, 10)
        L2 = deepcopy(L1)

        L2 .^= NF(2)
        @test L1 .^ NF(2) ≈ L2
    end
end 


@testset "LowerTriangularArray: GPU (JLArrays)" begin 
    # TODO: so far very basic GPU test, might integrate them into the other tests, as I already did with the broadcast test, but there are some key differences to avoid scalar indexing
    NF = Float32
    idims = (5,)

    L_cpu = randn(LowerTriangularArray{NF}, 10, 10, 5)

    # constructors/adapt
    L = adapt(JLArray, L_cpu)
    L2 = LowerTriangularArray(adapt(JLArray, L_cpu.data), 10, 10)
    @test all(L .== L2) 

    # getindex 
    @test typeof(L[1,:]) <: JLArray 
    for lm in SpeedyWeather.eachharmonic(L)
        @test Array(L[lm,:]) == L_cpu[lm,:]  
    end 

    # setindex! 
    A_test = JLArray(rand(NF,size(L_cpu,3)))
    L[1,:] = A_test
    @test L[1,:] == A_test

    # fill 
    fill!(L, 2)
    for lm in SpeedyWeather.eachharmonic(L2)
        @test all(L[lm, [Colon() for i=1:length(idims)]...] .== 2)
    end 

    # copy 
    L2 = deepcopy(L)
    @test all(L2 .== L)

    rand_array = JLArray(rand(NF,5))
    L2[1,:] = rand_array
    @test all(L[1,:] .== 2)     # should be a deep copy
    @test all(L2[1,:] .== rand_array)

    # rand + convert
    L3 = adapt(JLArray, randn(LowerTriangularArray{NF}, 10, 10, 5))
    L4 = convert(LowerTriangularArray{Float16,3,JLArray{Float16}}, L3)

    for lm in SpeedyWeather.eachharmonic(L, L3)
        @test all(Float16.(L3[lm, :]) .== L4[lm, :])
    end 
    
    # * 
    @test all((L+L) .== L*2)
    @test all((L-L) .== zero(L))

    # similar 
    @test size(similar(L)) == size(L)
    @test eltype(L) == eltype(similar(L, eltype(L)))

    @test (5, 7, 5) == size(similar(L, 5, 7, idims...))
    @test (5, 7, 5) == size(similar(L, (5, 7,  idims...)))
    @test similar(L) isa LowerTriangularArray

    # copyto! same size 
    L1 = randn(LowerTriangularArray{NF}, 10, 10, idims...)
    L1 = adapt(JLArray, L1)
    L2 = similar(L1)
    copyto!(L2, L1)

    @test all(L2 .== L1)
    
    # test the truncating copyto! function 
    # we can't do this with JLArrays, as they don't support mixed indexing with BitArrays
    # like Array and CuArray do
    # So, we do this with regular Array but with the _copyto_core! function that implements 
    # the core of this copyto! in a GPU compatible way, and is called by copyto! with CuArrays

    L1 = zeros(LowerTriangularArray{NF}, 33, 32, idims...);
    L2 = randn(LowerTriangularArray{NF}, 65, 64, idims...);
    L2T = spectral_truncation(L2,(size(L1)[1:2] .- 1)...)
    L3 = zeros(LowerTriangularArray{NF}, 33, 32, idims...);

    SpeedyWeather.LowerTriangularMatrices._copyto_core!(L1, L2, 1:33, 1:32)     # size of smaller matrix
    @test L1 == L2T

    # test that GPU and CPU method yield the same
    SpeedyWeather.LowerTriangularMatrices.copyto!(L3, L2, 1:33, 1:32)     # size of smaller matrix
    @test L1 == L3 

    SpeedyWeather.LowerTriangularMatrices._copyto_core!(L1, L2, 1:65, 1:64)     # size of bigger matrix
    @test L1 == L2T

    SpeedyWeather.LowerTriangularMatrices.copyto!(L3, L2, 1:65, 1:64)     # size of bigger matrix
    @test L1 == L3 

    SpeedyWeather.LowerTriangularMatrices._copyto_core!(L1, L2, 1:50, 1:50)     # in between
    @test L1 == L2T

    SpeedyWeather.LowerTriangularMatrices.copyto!(L3, L2, 1:50, 1:50)     # in between
    @test L3 == L1
end 

@testset "LowerTriangularArray: broadcast" begin 
    @testset for idims = ((), (5,), (5,5))
        @testset for NF in (Float16, Float32, Float64)
            @testset for ArrayType in (Array, JLArray)
                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = deepcopy(L1) 

                L2 .*= NF(5)
                @test all(L1 .* NF(5) .≈ L2) 

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = deepcopy(L1) 

                L2 ./= NF(5)
                @test all(L1 ./ NF(5) .≈ L2) 

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = deepcopy(L1)

                L2 .^= NF(2)
                @test all(L1 .^ NF(2) .≈ L2)

                # tests mirroring usage in dynamical core
                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L3 = deepcopy(L2)

                @. L2 += L1*NF(5)
                @test all(L2 .≈ (L3 + L1*NF(5)))

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))

                @. L3 = -L1 - L2
                @test all(L3 .≈ (-L1 - L2))

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L3 = deepcopy(L2)

                L2 .+= L1 
                L3 += L1
                @test all(L2 .≈ L3)

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = similar(L1)

                L2 .= NF(5) * L1
                @test all(L2 .≈ (L1*NF(5)))
            end
        end
    end
end 
