@testset "Grid indexing" begin
    for NF in (Float32,Float64)

        # with vector and resolution parameter provided
        L = FullClenshawGrid(randn(NF,96*47),24)        # L24 grid
        F = FullGaussianGrid(randn(NF,96*48),24)        # F24 grid
        O = OctahedralGaussianGrid(randn(NF,3168),24)   # O24 grid
        C = OctahedralClenshawGrid(randn(NF,3056),24)   # C24 grid
        H = HEALPixGrid(randn(NF,3072),16)              # H16 grid

        # without resolution parameter provided (inferred from vector length)
        L2 = FullClenshawGrid(randn(NF,96*47))          # L24 grid
        F2 = FullGaussianGrid(randn(NF,96*48))          # F24 grid
        O2 = OctahedralGaussianGrid(randn(NF,3168))     # O24 grid
        C2 = OctahedralClenshawGrid(randn(NF,3056))     # C24 grid
        H2 = HEALPixGrid(randn(NF,3072))                # H16 grid

        for (grid1,grid2) in zip([L,F,O,C,H],[L2,F2,O2,C2,H2])
            @test size(grid1) == size(grid2)
        end

        # getindex
        for ij in eachindex(L) L[ij] end
        for ij in eachindex(F) F[ij] end
        for ij in eachindex(O) O[ij] end
        for ij in eachindex(C) C[ij] end
        for ij in eachindex(H) H[ij] end

        # set index
        for ij in eachindex(L) L[ij] = 0 end
        for ij in eachindex(F) F[ij] = 0 end
        for ij in eachindex(O) O[ij] = 0 end
        for ij in eachindex(C) C[ij] = 0 end
        for ij in eachindex(H) H[ij] = 0 end

        @test all(L .== 0)
        @test all(F .== 0)
        @test all(O .== 0)
        @test all(C .== 0)
        @test all(H .== 0)
    end
end

@testset "Grid generators" begin
    for NF in (Float32,Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    )

            n = 32      # resolution parameter nlat_half/nside
            G1 = zeros(G,n)
            G2 = zero(G1)
            G3 = G(G2)

            @test G1 == G2
            @test G1 == G3
            @test all(G1 .== G2)
            @test all(G1 .== G3)

            # check that G2, G3 are deep copies
            G1[1] = 1
            @test G1[1] == 1
            @test G2[1] == 0
            @test G3[1] == 0
        end
    end
end

@testset "Grid indices" begin
    for G in (  FullClenshawGrid,
                FullGaussianGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                HEALPixGrid,
                )

        n = 32      # resolution parameter nlat_half/nside
        grid = zeros(G,n)

        for i in SpeedyWeather.eachring(grid)
            for ij in SpeedyWeather.each_index_in_ring(grid,i)
                grid[ij] += 1
            end
        end

        for ij in SpeedyWeather.eachgridpoint(grid)
            @test grid[ij] == 1
        end

        @test sum(grid) == SpeedyWeather.get_npoints(G,n)
    end
end
