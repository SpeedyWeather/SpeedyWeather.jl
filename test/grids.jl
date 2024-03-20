@testset "Grid indexing" begin
    for NF in (Float32, Float64)

        # with vector and resolution parameter provided
        L = FullClenshawGrid(randn(NF, 96*47), 24)        # L24 grid
        F = FullGaussianGrid(randn(NF, 96*48), 24)        # F24 grid
        O = OctahedralGaussianGrid(randn(NF, 3168), 24)   # O24 grid
        C = OctahedralClenshawGrid(randn(NF, 3056), 24)   # C24 grid
        H = HEALPixGrid(randn(NF, 3072), 32)              # H32 grid
        J = OctaHEALPixGrid(randn(NF, 4096), 32)             # J32 grid
        K = FullOctaHEALPixGrid(randn(NF, 128*63), 32)       # K32 grid

        # without resolution parameter provided (inferred from vector length)
        L2 = FullClenshawGrid(randn(NF, 96*47))          # L24 grid
        F2 = FullGaussianGrid(randn(NF, 96*48))          # F24 grid
        O2 = OctahedralGaussianGrid(randn(NF, 3168))     # O24 grid
        C2 = OctahedralClenshawGrid(randn(NF, 3056))     # C24 grid
        H2 = HEALPixGrid(randn(NF, 3072))                # H32 grid
        J2 = OctaHEALPixGrid(randn(NF, 4096))               # J32 grid
        K2 = FullOctaHEALPixGrid(randn(NF, 128*63))         # K32 grid

        for (grid1, grid2) in zip([L, F, O, C, H, J, K], [L2, F2, O2, C2, H2, J2, K2])
            @test size(grid1) == size(grid2)
        end

        # getindex
        for ij in eachindex(L) L[ij] end
        for ij in eachindex(F) F[ij] end
        for ij in eachindex(O) O[ij] end
        for ij in eachindex(C) C[ij] end
        for ij in eachindex(H) H[ij] end
        for ij in eachindex(J) J[ij] end
        for ij in eachindex(K) K[ij] end

        # set index
        for ij in eachindex(L) L[ij] = 0 end
        for ij in eachindex(F) F[ij] = 0 end
        for ij in eachindex(O) O[ij] = 0 end
        for ij in eachindex(C) C[ij] = 0 end
        for ij in eachindex(H) H[ij] = 0 end
        for ij in eachindex(J) J[ij] = 0 end
        for ij in eachindex(K) K[ij] = 0 end

        @test all(L .== 0)
        @test all(F .== 0)
        @test all(O .== 0)
        @test all(C .== 0)
        @test all(H .== 0)
        @test all(J .== 0)
        @test all(K .== 0)
    end
end

@testset "Grid generators" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            G1 = zeros(G{NF}, n)
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

            @test eltype(G1) == eltype(G2) == eltype(G3) == NF
        end
    end
end

@testset "Grid generators: ones" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            G1 = ones(G, n)
            @test all(G1 .== 1)
            @test eltype(G1) == Float64

            G2 = ones(G{NF}, n)
            @test all(G2 .== 1)
            @test eltype(G2) == NF
        end
    end
end

@testset "Grid generators: rand, randn" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            G1 = rand(G, n)
            @test eltype(G1) == Float64
            
            G1 = randn(G, n)
            @test eltype(G1) == Float64
            
            G1 = rand(G{NF}, n)
            @test eltype(G1) == NF
            
            G1 = randn(G{NF}, n)
            @test eltype(G1) == NF
        end
    end
end

@testset "Grid generators: undef" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            G1 = G(undef, n)
            @test eltype(G1) == Float64
            
            G1 = G{NF}(undef, n)
            @test eltype(G1) == NF
        end
    end
end

@testset "Grid indices" begin
    for G in (  FullClenshawGrid,
                FullGaussianGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                HEALPixGrid,
                OctaHEALPixGrid,
                FullHEALPixGrid,
                FullOctaHEALPixGrid,
                )

        n = 32      # resolution parameter nlat_half
        grid = zeros(G, n)

        # precompute indices and boundscheck
        rings = RingGrids.eachring(grid, grid)   

        for (j, ring) in enumerate(rings)
            for ij in ring
                grid[ij] += 1
            end
        end

        for ij in RingGrids.eachgridpoint(grid)
            @test grid[ij] == 1
        end

        @test sum(grid) == RingGrids.get_npoints(G, n)
    end
end

@testset "Ring indices" begin
    @testset for G in (  FullClenshawGrid,
                FullGaussianGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                HEALPixGrid,
                OctaHEALPixGrid,
                FullHEALPixGrid,
                FullOctaHEALPixGrid,
                )

        @testset for n in [8, 16, 24, 32]      # resolution parameter nlat_half
            grid = zeros(G, n)

            # precompute indices and boundscheck
            rings = SpeedyWeather.eachring(grid)   
            rings2 = [SpeedyWeather.each_index_in_ring(grid, j) for j in 1:SpeedyWeather.get_nlat(grid)]

            @test rings == rings2
        end
    end
end