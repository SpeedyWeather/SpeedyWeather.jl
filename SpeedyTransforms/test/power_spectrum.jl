@testset "Power spectrum" begin
    @testset for NF in (Float16, Float32, Float64)
        @testset for trunc in (31, 42)

            nlayers = 2

            spectrum = Spectrum(trunc)
            L = randn(complex(NF), spectrum, nlayers)

            # one more degree that shouldbe ignored
            spectrum2 = Spectrum(trunc, one_degree_more = true)
            L2 = randn(complex(NF), spectrum2, nlayers)

            # make trunc x trunc identical random numbers
            copyto!(L2, L)

            p = SpeedyTransforms.power_spectrum(L)
            p2 = SpeedyTransforms.power_spectrum(L2)

            @test p == p2
            @test eltype(p) == NF
            @test eltype(p2) == NF

            for k in 1:nlayers
                @test p[:, k] == SpeedyTransforms.power_spectrum(L[:, k])
            end

            # normalization should change things!
            @test p != SpeedyTransforms.power_spectrum(L, normalize = false)
        end
    end
end
