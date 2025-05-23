using JLArrays
using Adapt

@testset "Grid types" begin
    for G in (  
        FullClenshawGrid,
        FullGaussianGrid,
        OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        OctaminimalGaussianGrid,
        HEALPixGrid,
        OctaHEALPixGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid
        )

        full = RingGrids.full_grid_type(G)
        @test RingGrids.isfull(full)
        @test ~RingGrids.isreduced(full)
        @test full == RingGrids.full_grid_type(G)
        @test full <: RingGrids.AbstractFullGrid
        if ~(G <: RingGrids.AbstractFullGrid)
            @test G <: RingGrids.AbstractReducedGrid
        end
    end
end

@testset "Field types" begin
    # full ones
    for F in (  
        FullClenshawField,
        FullGaussianField,
        FullHEALPixField,
        FullOctaHEALPixField
        )

        @test RingGrids.isfull(F)
        @test ~RingGrids.isreduced(F)
        G = RingGrids.grid_type(F)
        @test RingGrids.isfull(G)
        @test RingGrids.field_type(RingGrids.grid_type(F)) <: F
    end

    for F in (
        OctahedralGaussianField,
        OctahedralClenshawField,
        OctaminimalGaussianField,
        HEALPixField,
        OctaHEALPixField,
    )

        @test RingGrids.isreduced(F)
        @test ~RingGrids.isfull(F)
        G = RingGrids.grid_type(F)
        @test RingGrids.isreduced(G)
        G = RingGrids.full_grid_type(F)
        @test RingGrids.isfull(G)
        @test RingGrids.field_type(RingGrids.grid_type(F)) <: F
    end
end

@testset "Grid indexing" begin
    for NF in (Float32, Float64)
        for Grid in (
            FullClenshawGrid,
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
            FullHEALPixGrid,
            FullOctaHEALPixGrid
            )
            for nlat_half in (4, 8, 16, 24, 32)
                
                npoints = RingGrids.get_npoints(Grid, nlat_half)
                grid = Grid(nlat_half)

                field1 = Field(grid)
                field2 = Field(zeros(npoints), grid)
                field3 = zeros(Grid{NF}, nlat_half)

                @test size(field1) == size(field2) == size(field3)
                @test_throws DimensionMismatch Field(zeros(npoints-1), grid)

                # getindex
                for ij in eachindex(field1) field1[ij] end

                # setindex
                for ij in eachindex(field1) field1[ij] = 0 end

                @test all(field1 .== 0)
                @test field1 == field3
            end
        end
    end
end

@testset "Field generators" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            grid = G(n)

            F1 = zeros(grid)
            F2 = zeros(grid, 1)
            F3 = zeros(NF, grid)
            F4 = zeros(NF, grid, 1)
            F5 = Field(zeros(NF, RingGrids.get_npoints(grid)), grid)
            F6 = Field(grid)
            F7 = Field(NF, grid)
            F8 = zero(F1)

            @test F1 == F2[:, 1]
            @test F1 == F3
            @test F1 == F4[:, 1]
            @test F1 == F5
            @test F1 == F6
            @test F1 == F7
            @test F1 == F8

            @test_broken all(F1 .== F2)
            @test all(F1 .== F3)
            @test_broken all(F1 .== F4)
            @test all(F1 .== F5)
            @test all(F1 .== F6)
            @test all(F1 .== F7)
            @test all(F1 .== F8)

            # check that they are deep copies
            F1[1] = 1
            @test F1[1] == 1
            @test F2[1] == 0
            @test F3[1] == 0
            @test F4[1] == 0
            @test F5[1] == 0
            @test F6[1] == 0
            @test F7[1] == 0
            @test F8[1] == 0

            @test eltype(F1) == eltype(F2) == eltype(F6) == eltype(F8)
            @test eltype(F3) == eltype(F4) == eltype(F5) == eltype(F7) == NF
        end
    end
end

@testset "Grid generators: ones" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
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

@testset "Field generators: rand, randn" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
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
            
            G1 = rand(NF, G, n)
            @test eltype(G1) == NF
            
            G1 = randn(NF, G, n)
            @test eltype(G1) == NF
        end
    end
end

@testset "Field generators: undef" begin
    for NF in (Float32, Float64)
        for G in (  FullClenshawGrid,
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    OctaminimalGaussianGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid,
                    FullHEALPixGrid,
                    FullOctaHEALPixGrid
                    )

            n = 4      # resolution parameter nlat_half
            F = RingGrids.field_type(G)
            field1 = F(undef, n)
            @test eltype(field1) == Float64
            
            field2 = F{NF}(undef, n)
            @test eltype(field2) == NF
        end
    end
end

@testset "FullGrids conversions to/from Arrays" begin 
    for idims in ((), (5,), (5,5))
        NF = Float64
        N = length(idims)+1
        data = rand(8,4, idims...)
        field = FullGaussianField(data, input_as=Matrix)
        @test Array(field, as=Matrix) == data

        data = rand(8,3, idims...)
        field = FullClenshawField(data, input_as=Matrix)
        @test Array(field, as=Matrix) == data

        data = rand(8,3, idims...)
        field = FullHEALPixField(data, input_as=Matrix)
        @test Array(field, as=Matrix) == data

        data = rand(8,3, idims...)
        field = FullOctaHEALPixField(data, input_as=Matrix)
        @test Array(field, as=Matrix) == data
    end   
end 

@testset "Grid indices" begin
    for G in (  FullClenshawGrid,
                FullGaussianGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                OctaminimalGaussianGrid,
                HEALPixGrid,
                OctaHEALPixGrid,
                FullHEALPixGrid,
                FullOctaHEALPixGrid,
                )

        n = 4      # resolution parameter nlat_half
        field = zeros(G, n)

        # precompute indices and boundscheck
        rings = RingGrids.eachring(field, field)   

        for (j, ring) in enumerate(rings)
            for ij in ring
                field[ij] += 1
            end
        end

        for ij in RingGrids.eachgridpoint(field)
            @test field[ij] == 1
        end

        @test sum(field) == RingGrids.get_npoints(G, n)
    end
end

@testset "Ring indices" begin
    @testset for G in (  FullClenshawGrid,
                FullGaussianGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                OctaminimalGaussianGrid,
                HEALPixGrid,
                OctaHEALPixGrid,
                FullHEALPixGrid,
                FullOctaHEALPixGrid,
                )

        @testset for n in [8, 16, 24, 32]      # resolution parameter nlat_half
            grid = G(n)

            # precompute indices and boundscheck
            rings = RingGrids.eachring(grid)   
            rings2 = [RingGrids.each_index_in_ring(grid, j) for j in 1:RingGrids.get_nlat(grid)]

            @test rings == rings2
        end
    end
end

@testset "Ring indices from fields" begin

    f1 = zeros(OctahedralGaussianGrid, 2)
    f2 = zeros(OctahedralGaussianGrid, 2, 1)    # matches above
    f3 = zeros(OctahedralGaussianGrid, 2, 2)    # matches horizontally only
    f4 = zeros(OctahedralClenshawGrid, 2)       # does not match above

    @test eachring(f1) == eachring(f1, f2) == eachring(f1, f2, f2, f1)
    @test eachring(f1) == eachring(f2, f3)
    @test_throws DimensionMismatch eachring(f1, f4)
    @test_throws DimensionMismatch eachring(f2, f4)
    @test_throws DimensionMismatch eachring(f3, f4)

    @test RingGrids.fields_match(f1, f3) == false
    @test RingGrids.fields_match(f2, f3) == false
    @test RingGrids.fields_match(f1, f3, horizontal_only=true)
    @test RingGrids.fields_match(f2, f3, horizontal_only=true)
end

@testset "Ring indices from grids" begin
    g1 = OctahedralGaussianGrid(2)
    g2 = OctahedralClenshawGrid(2)

    @test eachring(g1) != eachring(g2)
    @test eachring(g1) == eachring(g1, g1)
    @test_throws DimensionMismatch eachring(g1, g2)
    @test_throws DimensionMismatch eachring(g2, g1)

    @test RingGrids.grids_match(g1, g1)
    @test RingGrids.grids_match(g1, g2) == false
end

@testset "Field broadcasting" begin
    n = 2
    @testset for F in ( FullClenshawField,
                        FullGaussianField,            # don't test all to speed up CI
                        OctahedralGaussianField,
                        # OctahedralClenshawGrid,
                        # OctaminimalGaussianGrid,
                        HEALPixField,
                        # OctaHEALPixGrid,
                        # FullHEALPixGrid,
                        # FullOctaHEALPixGrid,
                        )

        @test zeros(F, n) .+ 1 == ones(F, n)
        @test ones(F, n)  .- 1 == zeros(F, n)
        @test ones(F, n)/1 == ones(F, n)
        @test zeros(F, n) + ones(F, n) == ones(F, n)
        @test 2ones(F, n) == ones(F, n) + ones(F, n)

        # don't promote to Array
        for s in ((n,), (n, n), (n, n, n), (n, n, n, n))
            field = zeros(F, s...)
            @test (field + field) isa Field
            @test (field .+ field) isa Field
            @test (field - field) isa Field
            @test (field .- field) isa Field
            @test (field .* field) isa Field
            @test (field ./ field) isa Field
            @test 2field isa Field
        end

        # promote types, Field{Float16} -> Field{Float64} etc
        @test all(ones(F{Float16}, n)*2.0 .=== 2.0)
        @test all(ones(F{Float16}, n)*2f0 .=== 2f0)
        @test all(ones(F{Float32}, n)*2.0 .=== 2.0)

        # promote types across grids
        @test all(ones(F{Float16}, n) + ones(F{Float32}, n) .=== 2f0)
        # @test all(ones(G{Float16}, n) + ones(G{Float64}, n) .=== 2.0)
        # @test all(ones(G{Float32}, n) + ones(G{Float64}, n) .=== 2.0)

        # # promote types across grids
        # @test all(ones(G{Float16}, n) - ones(G{Float32}, n) .=== 0f0)
        # @test all(ones(G{Float16}, n) - ones(G{Float64}, n) .=== 0.0)
        # @test all(ones(G{Float32}, n) - ones(G{Float64}, n) .=== 0.0)

        # # promote types across grids
        # @test all(ones(G{Float16}, n) .* ones(G{Float32}, n) .=== 1f0)
        # @test all(ones(G{Float16}, n) .* ones(G{Float64}, n) .=== 1.0)
        # @test all(ones(G{Float32}, n) .* ones(G{Float64}, n) .=== 1.0)

        # # promote types across grids
        # @test all(ones(G{Float16}, n) ./ ones(G{Float32}, n) .=== 1f0)
        # @test all(ones(G{Float16}, n) ./ ones(G{Float64}, n) .=== 1.0)
        # @test all(ones(G{Float32}, n) ./ ones(G{Float64}, n) .=== 1.0)

        # # dimension mismatches, but still broadcasting
        # g1 = rand(G, n)
        # g2 = rand(G, n, 1)
        # g2_2 = rand(G, n, 2)
        # g3 = rand(G, n, 1, 1)

        # @test (g1 .* g2)[:,1] == (g1 .* g2[:,1])
        # @test (g2 .* g1)[:,1] == (g1 .* g2[:,1])
        # @test (g2 .* g3)[:,1,1] == (g2[:,1] .* g3[:,1,1])
        # @test (g1 .* g2_2)[:,1] == (g1 .* g2_2[:,1])
        # @test (g1 .* g2_2)[:,2] == (g1 .* g2_2[:,2])
    end 
end

@testset "N-dimensional indexing" begin
    m, n, p = 2, 3, 4
    @testset for G in ( FullClenshawGrid,
                        FullGaussianGrid,
                        OctahedralGaussianGrid,
                        OctahedralClenshawGrid,
                        OctaminimalGaussianGrid,
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

        @test SpeedyWeather.RingGrids.nonparametric_type(typeof(grid[:,1:2,1:2])) <: RingGrids.nonparametric_type(G)
        @test grid[:, 1:2, 1:2].data == grid.data[:, 1:2, 1:2]

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
                        OctaminimalGaussianArray,
                        HEALPixArray,
                        OctaHEALPixArray,
                        FullHEALPixArray,
                        FullOctaHEALPixArray,
                        )

        for s in ((n,), (n, n), (n, n, n), (n, n, n, n))
            grid = zeros(G, s...)

            for k in eachlayer(grid)
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

# needed when extension is not loaded (manual testing)
# RingGrids.nonparametric_type(::Type{<:JLArray}) = JLArray

@testset "AbstractField: GPU (JLArrays)" begin 
    NF = Float32
    @testset for Grid in ( 
        FullClenshawArray,
        # FullGaussianArray,            # don't test all for CI speedup
        OctahedralGaussianArray,
        # OctahedralClenshawArray,
        # OctaminimalGaussianArray,
        HEALPixArray,
        # OctaHEALPixArray,
        # FullHEALPixArray,
        # FullOctaHEALPixArray,
    )
        s = (2, 3, 4)
        ndims = length(s)

        G_cpu = randn(Grid{NF}, s...)

        # constructors/adapt
        G = adapt(JLArray, G_cpu)
        G2 = Grid(adapt(JLArray, G_cpu.data))
        @test G == G2

        # broadcasting doesn't escape
        @test G  + G isa Grid{NF, ndims, JLArray{NF, ndims}}
        @test G .+ G isa Grid{NF, ndims, JLArray{NF, ndims}}
        @test G_cpu  + G_cpu isa Grid{NF, ndims, Array{NF, ndims}}
        @test G_cpu .+ G_cpu isa Grid{NF, ndims, Array{NF, ndims}}

        # getindex 
        @test G[1, :, :] isa JLArray{NF, 2}
        @test G[:, 1, 1] isa Grid{NF, 1, JLArray{NF, 1}}
        for k in eachlayer(G)
            for (j, ring) in enumerate(eachring(G))
                @test G[ring, k] == adapt(JLArray, G_cpu[ring, k])
                @test Array(G[ring, k]) == G_cpu[ring, k]
            end
        end

        # setindex! 
        G_test = JLArray(rand(NF,s[3]))
        G[1, 1, :] .= G_test                # with .
        @test G[1, 1, :] == G_test

        G_test = JLArray(rand(NF,s[3]))
        G[1, 1, :] = G_test                 # without .
        @test G[1, 1, :] == G_test

        # with other grid {Array}
        G_test = rand(Grid, s[1])
        G[:, 1, 1] = G_test                 # conversion to float64 -> float32
        @test Array(G[:, 1, 1].data) ≈ G_test

        # with other grid {JLArray}
        G_test = adapt(JLArray,rand(Grid, s[1]))
        G[:, 1, 1] = G_test
        @test G[:, 1, 1].data ≈ G_test.data

        # fill 
        fill!(G, 2)
        @test all(G .== 2)
    end
end

@testset "Zonal mean" begin
    @testset for NF in (Int32, Int64, Float16, Float32, Float64)
        @testset for Grid in ( 
            FullClenshawArray,
            FullGaussianArray,
            OctahedralGaussianArray,
            OctahedralClenshawArray,
            OctaminimalGaussianArray,
            HEALPixArray,
            OctaHEALPixArray,
            FullHEALPixArray,
            FullOctaHEALPixArray,
        )
            nlat_half = 4
            npoints = RingGrids.get_npoints(Grid, nlat_half)
            grid = Grid{NF}(1:npoints, nlat_half)
        
            zm = zonal_mean(grid)

            for (j, m) in enumerate(zonal_mean(grid))
                @test m == sum(RingGrids.eachring(grid)[j]) / RingGrids.get_nlon_per_ring(grid, j)
            end
        end
    end
end