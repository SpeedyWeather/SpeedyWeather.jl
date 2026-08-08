@testset "Create transform at uncommon resolutions" begin
    @testset for truncation in rand(15:200, 10)
        spectrum = Spectrum(truncation)
        grid = FullClenshawGrid(SpeedyTransforms.get_nlat_half(truncation))
        S = SpectralTransform(spectrum, grid)
    end
end

@testset "roundup nlon for FFT" begin
    for i in 1:10
        @test 2^i == SpeedyTransforms.roundup_fft(2^i)
        @test 2^i * 3 == SpeedyTransforms.roundup_fft(2^i * 3)
        @test 2^i * 5 == SpeedyTransforms.roundup_fft(2^i * 5)
    end
    for n in 1:10
        i = rand(2:1000)
        @test i <= SpeedyTransforms.roundup_fft(i)
    end
end
