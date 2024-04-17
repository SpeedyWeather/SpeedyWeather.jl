using JLArrays
using Adapt

@testset "Grid indexing" begin
    for NF in (Float32, Float64)

        # with vector and resolution parameter provided
        L = FullClenshawGrid(randn(NF, 96*47), 24)          # L24 grid
        F = FullGaussianGrid(randn(NF, 96*48), 24)          # F24 grid
        O = OctahedralGaussianGrid(randn(NF, 3168), 24)     # O24 grid
        C = OctahedralClenshawGrid(randn(NF, 3056), 24)     # C24 grid
        H = HEALPixGrid(randn(NF, 3072), 32)                # H32 grid
        J = OctaHEALPixGrid(randn(NF, 4096), 32)            # J32 grid
        K = FullOctaHEALPixGrid(randn(NF, 128*63), 32)      # K32 grid

        # without resolution parameter provided (inferred from vector length)
        L2 = FullClenshawGrid(randn(NF, 96*47))             # L24 grid
        F2 = FullGaussianGrid(randn(NF, 96*48))             # F24 grid
        O2 = OctahedralGaussianGrid(randn(NF, 3168))        # O24 grid
        C2 = OctahedralClenshawGrid(randn(NF, 3056))        # C24 grid
        H2 = HEALPixGrid(randn(NF, 3072))                   # H32 grid
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
            rings = RingGrids.eachring(grid)   
            rings2 = [RingGrids.each_index_in_ring(grid, j) for j in 1:RingGrids.get_nlat(grid)]

            @test rings == rings2
        end
    end
end

@testset "Grid broadcasting" begin
    n = 2
    @testset for G in ( FullClenshawArray,
                        FullGaussianArray,
                        OctahedralGaussianArray,
                        OctahedralClenshawArray,
                        HEALPixArray,
                        OctaHEALPixArray,
                        FullHEALPixArray,
                        FullOctaHEALPixArray,
                        )

        @test zeros(G, n) .+ 1 == ones(G, n)
        @test ones(G, n)  .- 1 == zeros(G, n)
        @test ones(G, n)/1 == ones(G, n)
        @test zeros(G, n) + ones(G, n) == ones(G, n)
        @test 2ones(G, n) == ones(G, n) + ones(G, n)

        # don't promote to Array
        for s in ((n,), (n, n), (n, n, n), (n, n, n, n))
            grid = zeros(G, s...)
            @test (grid + grid) isa G
            @test (grid - grid) isa G
            @test (grid .* grid) isa G
            @test (grid ./ grid) isa G
            @test 2grid isa G
        end

        # promote types, Grid{Float16} -> Grid{Float64} etc
        @test all(ones(G{Float16}, n)*2.0 .=== 2.0)
        @test all(ones(G{Float16}, n)*2f0 .=== 2f0)
        @test all(ones(G{Float32}, n)*2.0 .=== 2.0)

        # promote types across grids
        @test all(ones(G{Float16}, n) + ones(G{Float32}, n) .=== 2f0)
        @test all(ones(G{Float16}, n) + ones(G{Float64}, n) .=== 2.0)
        @test all(ones(G{Float32}, n) + ones(G{Float64}, n) .=== 2.0)

        # promote types across grids
        @test all(ones(G{Float16}, n) - ones(G{Float32}, n) .=== 0f0)
        @test all(ones(G{Float16}, n) - ones(G{Float64}, n) .=== 0.0)
        @test all(ones(G{Float32}, n) - ones(G{Float64}, n) .=== 0.0)

        # promote types across grids
        @test all(ones(G{Float16}, n) .* ones(G{Float32}, n) .=== 1f0)
        @test all(ones(G{Float16}, n) .* ones(G{Float64}, n) .=== 1.0)
        @test all(ones(G{Float32}, n) .* ones(G{Float64}, n) .=== 1.0)

        # promote types across grids
        @test all(ones(G{Float16}, n) ./ ones(G{Float32}, n) .=== 1f0)
        @test all(ones(G{Float16}, n) ./ ones(G{Float64}, n) .=== 1.0)
        @test all(ones(G{Float32}, n) ./ ones(G{Float64}, n) .=== 1.0)
    end 
end

@testset "N-dimensional indexing" begin
    m, n, p = 2, 3, 4
    @testset for G in ( FullClenshawGrid,
                        FullGaussianGrid,
                        OctahedralGaussianGrid,
                        OctahedralClenshawGrid,
                        HEALPixGrid,
                        OctaHEALPixGrid,
                        FullHEALPixGrid,
                        FullOctaHEALPixGrid,
                        )

        grid = rand(G, m, n, p)
        @test grid[:, 1, 1] isa G
        @test grid[1] == grid.data[1]
        @test grid[1, 1, 1] == grid.data[1, 1, 1]

        @test grid[1:2, 1:2, 1:2] == grid.data[1:2, 1:2, 1:2]
        @test grid[1, 1, :] == grid.data[1, 1, :]

        idx = CartesianIndex((1, 2, 3))
        @test grid[idx] == grid.data[idx]
        
        ids = CartesianIndices((m, n, p))
        @test grid[ids] == grid.data[ids]
        @test grid[ids] isa Array
    end
end

@testset "Loop indexing" begin
    n = 2
    @testset for G in ( FullClenshawArray,
                        FullGaussianArray,
                        OctahedralGaussianArray,
                        OctahedralClenshawArray,
                        HEALPixArray,
                        OctaHEALPixArray,
                        FullHEALPixArray,
                        FullOctaHEALPixArray,
                        )

        for s in ((n,), (n, n), (n, n, n), (n, n, n, n))
            grid = zeros(G, s...)

            for k in eachgrid(grid)
                for (j, ring) in enumerate(eachring(grid))
                    for ij in ring
                        grid[ij, k] = 1
                    end
                end
            end
            @test all(grid .== 1)
        end
    end
end

@testset "AbstractGridArray: GPU (JLArrays)" begin 
    NF = Float32
    @testset for Grid in ( 
        FullClenshawArray,
        FullGaussianArray,
        OctahedralGaussianArray,
        OctahedralClenshawArray,
        HEALPixArray,
        OctaHEALPixArray,
        FullHEALPixArray,
        FullOctaHEALPixArray,
    )
        G_cpu = randn(Grid{NF}, 2, 3, 4)

        # constructors/adapt
        G = adapt(JLArray, G_cpu)
        G2 = Grid(adapt(JLArray, G_cpu.data))
        # @test G == G2
        # @test all(G .== G2) 

        # getindex 
        # @test G[1, :] isa JLArray
        # @test G[:, 1] isa Grid{NF, 1, JLArray{NF}}
        # for k in eachgrid(G)
        #     for (j, ring) in enumerate(eachring(G))
        #         @test Array(L[lm,:]) == L_cpu[lm,:]  
        #     end
        # end 

        # # setindex! 
        # A_test = JLArray(rand(NF,size(L_cpu,3)))
        # L[1,:] = A_test
        # @test L[1,:] == A_test

        # # fill 
        # fill!(L, 2)
        # for lm in SpeedyWeather.eachharmonic(L2)
        #     @test all(L[lm, [Colon() for i=1:length(idims)]...] .== 2)
        # end 
    end
end