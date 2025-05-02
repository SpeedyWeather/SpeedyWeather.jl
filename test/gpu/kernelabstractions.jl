using KernelAbstractions

@testset "KernelAbstractions tests" begin 

    # only on CPU (currently)
    arch = SpeedyWeather.CPU()

    @kernel function mul_test!(A, @Const(B), @Const(C))
        i, j = @index(Global, NTuple)
        A[i, j] = B[i, j] * C[i, j]
    end
    
    B = SpeedyWeather.DeviceArray(device_setup, rand(10, 10))
    C = SpeedyWeather.DeviceArray(device_setup, rand(10, 10))
    A = zero(B)

    SpeedyWeather.launch!(arch, typeof(A), mul_test!, size(B), A, B, C)

    @test A ≈ B .* C 
end 