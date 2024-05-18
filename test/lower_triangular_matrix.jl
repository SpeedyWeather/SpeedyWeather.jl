using JLArrays, Adapt
import Random
import SpeedyWeather.LowerTriangularMatrices: matrix_size

@testset "LowerTriangularMatrix" begin
    @testset for NF in (Float32, Float64)
        mmax = 32
        @testset for lmax = (mmax, mmax+1)
            A = randn(Complex{NF}, lmax, mmax)

            SpeedyTransforms.spectral_truncation!(A)

            L = LowerTriangularMatrix(A)

            @test matrix_size(L) == size(A)

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
                ind = LowerTriangularMatrices.lowertriangle_indices(A)
                ind = @. ~(ind)
                A[ind] .= zero(Complex{NF})

                L = LowerTriangularArray(A)

                @test matrix_size(L) == size(A)
                @test size(L)[2:end] == size(A)[3:end]
                @test size(L)[1] == SpeedyWeather.LowerTriangularMatrices.nonzeros(size(A,1), size(A,2))

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
    m, n = 5, 5
    A = randn(LowerTriangularMatrix, m, n)
    
    # @testset "getindex" begin
        @test_throws BoundsError A[m, n+1]  # outside of i, j range
        @test_throws BoundsError A[m+1, n]  # outside of i, j range

        mn = LowerTriangularMatrices.nonzeros(m, n)
        @test_throws BoundsError A[mn+1]    # outside of k range

        # with @inbounds accessing [1,2] should return [5]
        # because the j > 1 is skipped and ij2k(1, 2, 5) = 5
        # which isn't correct but would never be called without @inbounds
        f(A, i) = @inbounds A[i]             # wrap into function
        f(A, i, j) = @inbounds A[i, j]

        @test f(A, 1) == A[1]       # valid
        @test f(A, 2, 1) == A[2]    # valid
        # @test f(A, 1, 2) == A[m]    # invalid

        @test f(A, CartesianIndex(2, 1)) == A[2, 1]
        # @test f(A, CartesianIndex(1, 2)) == A[n, 1] # invalid, disable for CI (*)
    # end

    # @testset "setindex!" begin
        @test_throws BoundsError A[m+1, n] = 1  # invalid
        @test_throws BoundsError A[m, n+1] = 1  # invalid
        @test_throws BoundsError A[1, 2] = 1    # upper triangle

        mn = LowerTriangularMatrices.nonzeros(m, n)
        @test_throws BoundsError A[mn+1] = 1

        # with @inbounds accessing [1,2] should return [5]
        # because the j > 1 is skipped and ij2k(1, 2, 5) = 5
        # which isn't correct but would never be called without @inbounds
        g!(A, i) = @inbounds A[i] = 1               # wrap into function
        g!(A, i, j) = @inbounds A[i, j] = 1

        g!(A, 1)                                    # valid
        @test f(A, 1) == A[1] == 1

        g!(A, 2, 1)                                 # valid
        @test f(A, 2, 1) == A[2] == A[2,1] == 1

        # g!(A, 1, 2)                                 # invalid, disable for CI (*)
        # @test f(A, 1, 2) == A[n] == A[n,1] == 1
    # end

    # (*) github actions CI is by default set to _always_ checkbounds
    # regardless @inboounds annotations, so these tests were performed manually (and pass)
    # but are disabled to make github actions CI pass
end

@testset "LowerTriangularArray: @inbounds" begin
    m, n, p = 5, 5, 5
    A = randn(LowerTriangularArray, m, n, p)
    
    # @testset "getindex" begin
        @test_throws BoundsError A[m, n+1, p+1]  # outside or range
        @test_throws BoundsError A[m, n+1, p  ]  # outside or range
        @test_throws BoundsError A[m+1, n, p  ]  # outside or range

        mnp = LowerTriangularMatrices.nonzeros(m, n)*p
        @test_throws BoundsError A[mnp+1]       # outside of k range

        # with @inbounds accessing [1,2] should return [5]
        # because the j > 1 boundscheck is skipped and ij2k(1, 2, 5) = 5
        # which isn't correct but would never be called without @inbounds
        f(A, i) = @inbounds A[i]            # wrap into function
        f(A, i, j) = @inbounds A[i, j]
        f(A, i, j, k) = @inbounds A[i, j, k]

        @test f(A, 1) == A[1]               # valid
        @test f(A, 2, 1) == A[2]            # valid
        # @test f(A, 1, 2, 3) == A[m, 1, 3]   # invalid, disable for CI (*)
        @test f(A, 1, 1, 1) == A[1, 1, 1]   # valid

        @test f(A, CartesianIndex(3)) == A[3]
        @test f(A, CartesianIndex(2, 1)) == A[2, 1]
        @test f(A, CartesianIndex(2, 1, 1)) == A[2, 1]
        # @test f(A, CartesianIndex(1, 2, 1)) == A[n, 1]  # invalid, disable for CI (*)
    # end

    # @testset "setindex!" begin
        @test_throws BoundsError A[m+1, n, p  ] = 1  # invalid
        @test_throws BoundsError A[m, n+1, p  ] = 1  # invalid
        @test_throws BoundsError A[m, n+1, p+1] = 1  # invalid
        @test_throws BoundsError A[1, 2,   1  ] = 1  # upper triangle

        mnp = LowerTriangularMatrices.nonzeros(m, n)*p
        @test_throws BoundsError A[mnp+1] = 1

        # with @inbounds accessing [1,2] should return [5]
        # because the j > 1 is skipped and ij2k(1, 2, 5) = 5
        # which isn't correct but would never be called without @inbounds
        g!(A, i) = @inbounds A[i] = 1               # wrap into function
        g!(A, i, j) = @inbounds A[i, j] = 1
        g!(A, i, j, k) = @inbounds A[i, j, k] = 1

        g!(A, 1)                                    # valid
        @test f(A, 1) == A[1] == 1

        g!(A, 2, 1)                                 # valid
        @test f(A, 2, 1) == A[2] == A[2,1] == 1

        g!(A, 1, 2)                                 # valid
        @test f(A, 1, 2) == A[1, 2] == 1

        # g!(A, 1, 2, 1)                              # invalid, disable for CI (*)
        # @test f(A, 1, 2, 1) == A[n] == A[n, 1] == A[n, 1, 1] == 1
    # end
end

@testset "4D LowerTriangularArray: @inbounds" begin
    A = randn(LowerTriangularArray, 33, 32, 1, 1)
    
    @testset "getindex" begin
        @test_throws BoundsError A[34, 32, 1, 1]    # outside of i, j range
        @test_throws BoundsError A[561, 1, 1]       # outside of k range
        @test_throws BoundsError A[33, 32, 2, 1]    # outside of 3rd dim range
        @test_throws BoundsError A[33, 32, 1, 2]    # outside of 4th dim range
    end

    @testset "setindex!" begin
        @test_throws BoundsError A[34, 32, 1, 1] = 1
        @test_throws BoundsError A[561, 1, 1] = 1
        @test_throws BoundsError A[33, 32, 2, 1] = 1 
        @test_throws BoundsError A[33, 32, 1, 2] = 1 
    end
end

@testset "Zeros, ones, rand, and randn constructors" begin
    for f in (ones, zeros, rand, randn)
        s = (5, 5)
        
        # for 2D doesn't matter whether you say Matrix or Array, size is determined by s
        L = f(LowerTriangularMatrix, s...)
        L2 = f(LowerTriangularArray, s...)
        @test typeof(L) == typeof(L2)
        @test size(L) == size(L2)
        
        L = f(LowerTriangularMatrix{Float16}, s...)
        L2 = f(LowerTriangularArray{Float16}, s...)
        @test typeof(L) == typeof(L2)
        @test size(L) == size(L2)
        
        JL = adapt(JLArray, L)
        JL2 = adapt(JLArray, L2)
        @test typeof(JL) == typeof(JL2) == typeof(zero(JL))
        @test size(JL) == size(JL2) == size(zero(JL))
        
        L = f(LowerTriangularMatrix{Float16}, s...)
        L2 = f(LowerTriangularArray{Float16, 1, Vector{Float16}}, s...)
        JL = f(LowerTriangularArray{Float16, 1, JLArray{Float16, 1}}, s...)
        @test typeof(L) == typeof(L2)
        @test size(L) == size(L2)
        @test typeof(L) != typeof(JL)
        @test size(L) == size(L2)

        # higher dims
        for s in ((2, 3, 4), (2, 3, 4, 5))
            N = length(s)
            Random.seed!(123)
            L =  f(LowerTriangularArray{Float16, N-1,   Array{Float16, N-1}}, s...)
            Random.seed!(123)
            JL = f(LowerTriangularArray{Float16, N-1, JLArray{Float16, N-1}}, s...)
            JL2 = adapt(JLArray, L)
            @test all(JL2 .== JL)   # equality via broadcasting
            @test JL2 == JL         # checks for type and data equality
        end
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
                L3 = convert(LowerTriangularArray{Float16, 1+length(idims), Array{Float16,1+length(idims)}}, L)
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

        @test L+L isa LowerTriangularMatrix
        @test 2L isa LowerTriangularMatrix

        @test eachindex(L) == eachindex(L.data)
        @test eachindex(L, L) == eachindex(L.data, L.data)

        @test size(similar(L)) == size(L)
        @test eltype(L) == eltype(similar(L, eltype(L)))

        @test (5, 7) == matrix_size(similar(L, 5, 7))
        @test (5, 7) == matrix_size(similar(L, (5, 7)))
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

            @test L+L isa LowerTriangularArray
            @test 2L isa LowerTriangularArray

            @test eachindex(L) == eachindex(L.data)
            @test eachindex(L, L) == eachindex(L.data, L.data)

            @test size(similar(L)) == size(L)
            @test eltype(L) == eltype(similar(L, eltype(L)))

            @test (5, 7, idims...) == matrix_size(similar(L, 5, 7, idims...))
            @test (5, 7, idims...) == matrix_size(similar(L, (5, 7,  idims...)))
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
        L2T = spectral_truncation(L2,(matrix_size(L1) .- 1)...)

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

            @test L2 == L1
            
            # copyto! Array 
            M = zeros(NF, 10, 10, idims...)
            copyto!(M, L1)

            ind = SpeedyWeather.LowerTriangularMatrices.lowertriangle_indices(M)
            not_ind = @. ~(ind)

            @test all(M[not_ind] .== zero(NF))
            @test LowerTriangularArray(M) == L1

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
            L2T = spectral_truncation(L2,(matrix_size(L1)[1:2] .- 1)...)

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
    A_test = JLArray(rand(NF,size(L_cpu,2)))
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
    L4 = convert(LowerTriangularArray{Float16,2,JLArray{Float16,2}}, L3)

    for lm in SpeedyWeather.eachharmonic(L, L3)
        @test all(Float16.(L3[lm, :]) .== L4[lm, :])
    end 
    
    # * 
    @test all((L+L) .== L*2)
    @test all((L-L) .== zero(L))

    # similar 
    @test size(similar(L)) == size(L)
    @test eltype(L) == eltype(similar(L, eltype(L)))

    @test (5, 7, 5) == matrix_size(similar(L, 5, 7, idims...))
    @test (5, 7, 5) == matrix_size(similar(L, (5, 7,  idims...)))
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
    L2T = spectral_truncation(L2,(matrix_size(L1)[1:2] .- 1)...)
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

                L2 .*= 5
                @test 5L1 == L2

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = deepcopy(L1) 

                L2 ./= 5
                @test L1/5 == L2

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = deepcopy(L1)

                L2 .^= 2
                @test L1.^2 == L2

                # tests mirroring usage in dynamical core
                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L3 = deepcopy(L2)

                @. L2 += 5L1
                @test L2 == L3 + 5L1

                SpeedyWeather.LowerTriangularMatrices.add!(L3,5L1)
                @test L3 == L2 

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))

                @. L3 = -L1 - L2
                @test L3 == -L1 - L2

                L4 = adapt(ArrayType, zeros(LowerTriangularArray{NF}, 10, 10, idims...))
                SpeedyWeather.LowerTriangularMatrices.add!(L4,-L1,-L2)
                @test L4 == L3

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L3 = deepcopy(L2)

                L2 .+= L1 
                L3 += L1
                @test L2 == L3

                L1 = adapt(ArrayType, randn(LowerTriangularArray{NF}, 10, 10, idims...))
                L2 = similar(L1)

                L2 .= 5L1
                @test L2 == 5L1

                # add! test here
            end
        end
    end
end 
