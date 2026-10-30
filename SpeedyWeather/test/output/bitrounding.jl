@testset "Bitrounding preserves NaN missing values" begin
    for (T, U, all_mantissa_bits_nan) in ((Float32, UInt32, 0x7fffffff), (Float64, UInt64, 0x7fffffffffffffff))
        nan = reinterpret(T, all_mantissa_bits_nan)    # NaN as produced on CUDA GPUs
        @test isnan(nan)

        # without canonicalization the rounding carry flips the sign bit, turning the NaN into -0.0
        A = [nan, T(1.2345), T(-6.789)]
        SpeedyWeather.round!(A, 7)
        @test !isnan(A[1])

        # with canonicalization, NaNs survive rounding and finite values are rounded as before
        A = [nan, T(1.2345), T(-6.789)]
        B = copy(A)
        SpeedyWeather.canonicalize_nans!(A)
        @test isnan(A[1])
        @test A[2:3] == B[2:3]
        SpeedyWeather.round!(A, 7)
        SpeedyWeather.round!(B, 7)
        @test isnan(A[1])
        @test A[2:3] == B[2:3]
    end
end
