using LowerTriangularArrays
using SpeedyWeatherInternals.ArrayDimensions: LM, LMZ, LMT, LMZT,
    DimensionsWithTime, DimensionsWithVertical, DimensionsWithTimeAndVertical

@testset "LowerTriangularArray: trait-based dispatch" begin
    arch = LowerTriangularArrays.CPU()
    spectrum = Spectrum(10, 10, architecture=arch)

    # Create test arrays with different dimensions
    data_2d = rand(ComplexF32, 55)  # nonzeros for T=10, T=10
    data_3d_vertical = rand(ComplexF32, 55, 5)  # + 1 vertical level
    data_3d_time = rand(ComplexF32, 55, 5)  # + 1 time level
    data_4d = rand(ComplexF32, 55, 5, 10)  # + vertical and time

    # 2D array (LM)
    lta_lm = LowerTriangularArray(data_2d, spectrum, LM())
    @test lta_lm isa LowerTriangularArray
    @test !(lta_lm isa LowerTriangularArrayWithTime)
    @test !(lta_lm isa LowerTriangularArrayWithVertical)

    # 3D with vertical (LMZ)
    lta_lmz = LowerTriangularArray(data_3d_vertical, spectrum, LMZ())
    @test lta_lmz isa LowerTriangularArray
    @test !(lta_lmz isa LowerTriangularArrayWithTime)
    @test lta_lmz isa LowerTriangularArrayWithVertical

    # 3D with time (LMT)
    lta_lmt = LowerTriangularArray(data_3d_time, spectrum, LMT())
    @test lta_lmt isa LowerTriangularArray
    @test lta_lmt isa LowerTriangularArrayWithTime
    @test !(lta_lmt isa LowerTriangularArrayWithVertical)

    # 4D with both time and vertical (LMZT)
    lta_lmzt = LowerTriangularArray(data_4d, spectrum, LMZT())
    @test lta_lmzt isa LowerTriangularArray
    @test lta_lmzt isa LowerTriangularArrayWithTime
    @test lta_lmzt isa LowerTriangularArrayWithVertical
    @test lta_lmzt isa LowerTriangularArrayWithTimeAndVertical
end
