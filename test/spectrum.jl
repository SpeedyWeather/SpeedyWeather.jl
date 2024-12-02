@testset "Power spectrum" begin
    @testset for T in (Float16, Float32, Float64)
        @testset for trunc in (31, 42)

            ks = 2

            L = randn(LowerTriangularMatrix{complex(T)}, trunc+1, trunc+1, ks)

            # one more degree that shouldbe ignored
            L2 = randn(LowerTriangularMatrix{complex(T)}, trunc+2, trunc+1, ks)

            # make trunc x trunc identical random numbers
            copyto!(L2, L)

            p = SpeedyTransforms.power_spectrum(L)
            p2 = SpeedyTransforms.power_spectrum(L2)

            @test p == p2
            @test eltype(p) == T
            @test eltype(p2) == T

            for k in 1:ks
                @test p[:, k] == SpeedyTransforms.power_spectrum(L[:, k])
            end

            # normalization should change things!
            @test p != SpeedyTransforms.power_spectrum(L, normalize=false)
        end
    end
end