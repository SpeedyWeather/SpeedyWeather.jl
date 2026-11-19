@testset "get_step on plain arrays (as inside GPU kernels)" begin
    # a vector has no step dimension: it is returned as it is (step 1), not as a 2-index view (fails on Metal)
    v = rand(Float32, 10)
    @test SpeedyWeather.get_step(v, 1) === v

    # matrices and 3-D arrays: the last dimension is the step
    m = rand(Float32, 10, 2)
    @test SpeedyWeather.get_step(m, 2) == m[:, 2]
    a = rand(Float32, 10, 3, 2)
    @test SpeedyWeather.get_step(a, 1) == a[:, :, 1]
end
