using RingGrids
using SpeedyWeatherInternals.ArrayDimensions: XY, XYZ, XYT, XYZT, ZXY, ZXYT, LM, LMZ, LMT, LMZT,
    Dimensions2D, Dimensions3D, Dimensions4D,
    DimensionsWithTime, DimensionsWithVertical, DimensionsWithTimeAndVertical

@testset "AbstractField: 2D/3D/4D dispatch" begin
    arch = CPU()
    grid = HEALPixGrid(6, arch)

    # 2D field
    field_2d = rand(Float32, grid)
    @test field_2d isa AbstractField2D

    # 3D field with vertical
    field_3d_vertical = Field(zeros(Float32, grid, 5), grid, XYZ())
    @test field_3d_vertical isa AbstractField3D

    # 3D field with time
    field_3d_time = Field(zeros(Float32, grid, 5), grid, XYT())
    @test field_3d_time isa AbstractField3D

    # 4D field
    field_4d = Field(zeros(Float32, grid, 5, 10), grid, XYZT())
    @test field_4d isa AbstractField4D
end

@testset "AbstractField: trait-based dispatch" begin
    arch = CPU()
    grid = HEALPixGrid(6, arch)

    # Fields with vertical dimension (XYZ, XYZT)
    field_xyz = Field(zeros(Float32, grid, 5), grid, XYZ())
    field_xyzt = Field(zeros(Float32, grid, 5, 10), grid, XYZT())

    @test field_xyz isa AbstractFieldWithVertical
    @test field_xyzt isa AbstractFieldWithVertical

    # Fields with time dimension (XYT, XYZT)
    field_xyt = Field(zeros(Float32, grid, 5), grid, XYT())
    field_xyzt2 = Field(zeros(Float32, grid, 5, 10), grid, XYZT())

    @test field_xyt isa AbstractFieldWithTime
    @test field_xyzt2 isa AbstractFieldWithTime

    # Fields with both time and vertical (XYZT)
    field_xyzt3 = Field(zeros(Float32, grid, 5, 10), grid, XYZT())

    @test field_xyzt3 isa AbstractFieldWithTimeAndVertical

    # Negative tests: fields without time should not match AbstractFieldWithTime
    field_xy = rand(Float32, grid)
    field_xyz_no_time = Field(zeros(Float32, grid, 5), grid, XYZ())

    @test !(field_xy isa AbstractFieldWithTime)
    @test !(field_xyz_no_time isa AbstractFieldWithTime)
end
