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



    @kernel ∇²_kernel!(∇²alms, alms, @Const(mode_func))
        m, k = @index(Global, NTuple)


        ∇²alms[lm, k] = mode_func(∇²alms[lm, k], alms[lm, k]*eigenvalues[l])
    end 
end 