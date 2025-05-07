using KernelAbstractions
import SpeedyWeather: on_architecture
@testset "KernelAbstractions tests" begin 

    # To-Do write tests for each type of dims_type in the kernel launching util

    NF= Float32
    alms = rand(LowerTriangularArray{Complex{NF}},33, 32, 2)
    alms2 = copy(alms)
    alms3 = copy(alms)

    S = SpectralTransform(alms)

    SpeedyWeather.SpeedyTransforms.∇²_KA!(alms2, alms, S);

    SpeedyWeather.SpeedyTransforms.∇²!(alms3, alms, S);

    @test alms3 ≈ alms2
end 