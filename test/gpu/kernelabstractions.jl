using KernelAbstractions
import SpeedyWeather: on_architecture
@testset "KernelAbstractions tests" begin 

    # only on CPU (currently)
    arch = SpeedyWeather.CPU()

    @kernel function mul_test!(A, @Const(B), @Const(C))
        i, j = @index(Global, NTuple)
        A[i, j] = B[i, j] * C[i, j]
    end
    
    B = on_architecture(arch, rand(10, 10))
    C = on_architecture(arch, rand(10, 10))
    A = zero(B)

    SpeedyWeather.launch!(arch, typeof(A), size(B), mul_test!, A, B, C)

    @test A ≈ B .* C 

    NF= Float32
    alms = rand(LowerTriangularArray{Complex{NF}},33, 32, 2)
    alms2 = copy(alms)
    alms3 = copy(alms)

    S = SpectralTransform(alms)

    SpeedyWeather.SpeedyTransforms.∇²_KA!(alms2, alms, S);

    SpeedyWeather.SpeedyTransforms.∇²!(alms3, alms, S);

    @test alms3 ≈ alms2
end 