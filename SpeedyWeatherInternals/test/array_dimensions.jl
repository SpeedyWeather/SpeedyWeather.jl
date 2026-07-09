using SpeedyWeatherInternals.ArrayDimensions: XY, XYZ, XYT, XYZT, ZXY, ZXYT, LM, LMZ, LMT, LMZT,
    hastime, hasvertical,
    Dimensions2D, Dimensions3D, Dimensions4D,
    DimensionsWithTime, DimensionsWithVertical, DimensionsWithTimeAndVertical

@testset "ArrayDimensions: ndims" begin
    @test ndims(XY()) == 1
    @test ndims(XYZ()) == 2
    @test ndims(XYT()) == 2
    @test ndims(XYZT()) == 3
    @test ndims(LM()) == 1
    @test ndims(LMZ()) == 2
    @test ndims(LMT()) == 2
    @test ndims(LMZT()) == 3
end

@testset "ArrayDimensions: hastime, hasvertical" begin
    @test !hastime(XY()) && !hasvertical(XY())
    @test !hastime(XYZ()) && hasvertical(XYZ())
    @test hastime(XYT()) && !hasvertical(XYT())
    @test hastime(XYZT()) && hasvertical(XYZT())
    @test !hastime(LM()) && !hasvertical(LM())
    @test !hastime(LMZ()) && hasvertical(LMZ())
    @test hastime(LMT()) && !hasvertical(LMT())
    @test hastime(LMZT()) && hasvertical(LMZT())
end

@testset "ArrayDimensions: getindex drops dims" begin
    # 3D → 2D: integer index on non-horizontal dim
    @test XYZ()[:, 1] isa XY
    @test XYT()[:, 1] isa XY
    @test LMZ()[:, 1] isa LM
    @test LMT()[:, 1] isa LM

    # 4D → 3D: integer on vertical keeps time, integer on time keeps vertical
    @test XYZT()[:, 1] isa XYT          # drop vertical
    @test XYZT()[:, 1:3, 1] isa XYZ     # drop time
    @test LMZT()[:, 1] isa LMT          # drop vertical
    @test LMZT()[:, 1:3, 1] isa LMZ     # drop time

    # 4D → 2D: integer on both
    @test XYZT()[:, 1, 1] isa XY
    @test LMZT()[:, 1, 1] isa LM

    # range index preserves dims
    @test XYZ()[:, 1:3] isa XYZ
    @test LMZ()[:, 1:3] isa LMZ
    @test XYZT()[:, 1:3] isa XYZT
    @test LMZT()[:, 1:3] isa LMZT
end

@testset "ArrayDimensions: 2D/3D/4D union dispatch" begin
    # 2D types
    @test XY() isa Dimensions2D
    @test LM() isa Dimensions2D

    # 3D types
    @test XYZ() isa Dimensions3D
    @test XYT() isa Dimensions3D
    @test ZXY() isa Dimensions3D
    @test LMZ() isa Dimensions3D
    @test LMT() isa Dimensions3D

    # 4D types
    @test XYZT() isa Dimensions4D
    @test ZXYT() isa Dimensions4D
    @test LMZT() isa Dimensions4D
end

@testset "ArrayDimensions: trait-based union dispatch" begin
    # DimensionsWithTime: XYT, XYZT, LMT, LMZT, ZXYT
    @test XYT() isa DimensionsWithTime
    @test XYZT() isa DimensionsWithTime
    @test LMT() isa DimensionsWithTime
    @test LMZT() isa DimensionsWithTime
    @test ZXYT() isa DimensionsWithTime

    # DimensionsWithVertical: XYZ, XYZT, LMZ, LMZT, ZXY, ZXYT
    @test XYZ() isa DimensionsWithVertical
    @test XYZT() isa DimensionsWithVertical
    @test LMZ() isa DimensionsWithVertical
    @test LMZT() isa DimensionsWithVertical
    @test ZXY() isa DimensionsWithVertical
    @test ZXYT() isa DimensionsWithVertical

    # DimensionsWithTimeAndVertical: XYZT, ZXYT, LMZT
    @test XYZT() isa DimensionsWithTimeAndVertical
    @test ZXYT() isa DimensionsWithTimeAndVertical
    @test LMZT() isa DimensionsWithTimeAndVertical

    # Negative tests: types not in the unions
    @test !(XY() isa DimensionsWithTime)
    @test !(XY() isa DimensionsWithVertical)
    @test !(XYZ() isa DimensionsWithTime)
    @test !(XYT() isa DimensionsWithVertical)
    @test !(LM() isa DimensionsWithTime)
    @test !(LM() isa DimensionsWithVertical)
end
