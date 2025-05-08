using KernelAbstractions
import SpeedyWeather: on_architecture
@testset "KernelAbstractions tests" begin 

    # To-Do write tests for each type of dims_type in the kernel launching util, 
    # the tests currently below will be removed when the KA becomes the only one
    using SpeedyWeather, BenchmarkTools, Test
    # ∇²!

    NF= Float32
    alms = rand(LowerTriangularArray{Complex{NF}},33, 32, 2)
    #alms = rand(LowerTriangularArray{Complex{NF}},33, 32)

    alms2 = copy(alms)
    alms3 = copy(alms)
    alms4 = copy(alms)

    S = SpectralTransform(alms)

    SpeedyWeather.SpeedyTransforms.∇²!(alms2, alms, S);
    SpeedyWeather.SpeedyTransforms.∇²_KA!(alms3, alms, S);

    @test alms3 ≈ alms2

    # Divergence

    alms = rand(LowerTriangularArray{Complex{NF}},33, 32, 2)
    alms2 = rand(LowerTriangularArray{Complex{NF}},33, 32, 2)

    alms3 = copy(alms)
    alms4 = copy(alms)

    SpeedyWeather.SpeedyTransforms.divergence!(alms3, alms, alms2, S)
    SpeedyWeather.SpeedyTransforms.divergence_KA!(alms4, alms, alms2, S)

    @test alms4 ≈ alms3

end 