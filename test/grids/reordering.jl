@testset "Reordering identity" begin

    @testset for nlat_half in (4, 8, 16, 32, 64)
        grid = OctaHEALPixGrid(nlat_half)
        N = RingGrids.get_npoints(grid)

        @test [RingGrids.rcq2ring(RingGrids.ring2rcq(ij, grid)..., grid) for ij in 1:N] == 1:N
        @test [RingGrids.nest2ring(RingGrids.ring2nest(ij, grid)..., grid) for ij in 1:N] == 1:N
        @test [RingGrids.rcq2nest(RingGrids.nest2rcq(ij, grid)..., grid) for ij in 1:N] == 1:N

        field = rand(grid)
        field2 = RingGrids.ring_order(RingGrids.nested_order(field))
        field3 = RingGrids.nested_order(RingGrids.ring_order(field))
        @test field == field2
        @test field == field3

        field4 = RingGrids.matrix_order(field)
        M = Matrix(field4)
        @test size(M) == (2nlat_half, 2nlat_half)
    end
end

@testset "Interleaving bits" begin
    @testset for T in (Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64)
        m = min(typemax(T) >> (sizeof(T)*4), 2^18 - 1)  # limit to 18 bits for performance
        for i in 0:m
            ui = T(i)
            @test RingGrids.deinterleave(RingGrids.interleave_with_zeros(ui)) == ui
        end

        # random bitpatterns too
        for _ in 1:1024
            ui = T(rand(0:m))
            @test RingGrids.deinterleave(RingGrids.interleave_with_zeros(ui)) == ui
        end
    end
end