using KernelAbstractions

@testset "KernelAbstractions tests" begin 

    # only on CPU (currently)
    device_setup = SpeedyWeather.DeviceSetup(SpeedyWeather.CPUDevice())

    @kernel function mul_test!(A, @Const(B), @Const(C))
        i, j = @index(Global, NTuple)
        A[i,j] = B[i,j] * C[i,j]
    end
    
    B = SpeedyWeather.DeviceArray(device_setup, rand(10,10))
    C = SpeedyWeather.DeviceArray(device_setup, rand(10,10))
    A = zero(B)

    ev = SpeedyWeather.launch_kernel!(device_setup, mul_test!, size(B), A, B, C)
    wait(ev)

    @test A ≈ B .* C 
end 