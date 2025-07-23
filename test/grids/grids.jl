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
                field3 = zeros(NF, Grid, nlat_half)

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

            @test all(F1 .== F2)
            @test all(F2 .== F1)
            @test all(F1 .== F3)
            @test all(F1 .== F4)
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

@testset "Field generators: ones" begin
    for NF in (Float32, Float64)
        for F in (  FullClenshawField,
                    FullGaussianField,
                    OctahedralGaussianField,
                    OctahedralClenshawField,
                    OctaminimalGaussianField,
                    HEALPixField,
                    OctaHEALPixField,
                    FullHEALPixField,
                    FullOctaHEALPixField,
                    )

            n = 4      # resolution parameter nlat_half
            f1 = ones(F, n)
            @test all(f1 .== 1)
            @test eltype(f1) == Float64

            f2 = ones(F{NF}, n)
            @test all(f2 .== 1)
            @test eltype(f2) == NF
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
        for F in (  FullClenshawField,
                    FullGaussianField,
                    OctahedralGaussianField,
                    OctahedralClenshawField,
                    OctaminimalGaussianField,
                    HEALPixField,
                    OctaHEALPixField,
                    FullHEALPixField,
                    FullOctaHEALPixField,
                    )

            n = 4      # resolution parameter nlat_half
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

        for n in [8, 16, 24, 32]      # resolution parameter nlat_half
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
        @test all(ones(F{Float16}, n) + ones(F{Float64}, n) .=== 2.0)
        @test all(ones(F{Float32}, n) + ones(F{Float64}, n) .=== 2.0)

        # promote types across grids
        @test all(ones(F{Float16}, n) - ones(F{Float32}, n) .=== 0f0)
        @test all(ones(F{Float16}, n) - ones(F{Float64}, n) .=== 0.0)
        @test all(ones(F{Float32}, n) - ones(F{Float64}, n) .=== 0.0)

        # promote types across grids
        # f3 = ones(F{Float16}, n) .* ones(F{Float32}, n)
        @test all((ones(F{Float16}, n) .* ones(F{Float32}, n)) .=== 1f0)
        @test all(ones(F{Float16}, n) .* ones(F{Float64}, n) .=== 1.0)
        @test all(ones(F{Float32}, n) .* ones(F{Float64}, n) .=== 1.0)

        # promote types across grids
        @test all(ones(F{Float16}, n) ./ ones(F{Float32}, n) .=== 1f0)
        @test all(ones(F{Float16}, n) ./ ones(F{Float64}, n) .=== 1.0)
        @test all(ones(F{Float32}, n) ./ ones(F{Float64}, n) .=== 1.0)

        # dimension mismatches, but still broadcasting
        field1 = rand(F, n)
        field2 = rand(F, n, 1)
        field2_2 = rand(F, n, 2)
        field3 = rand(F, n, 1, 1)

        @test (field1 .* field2)[:,1]   == (field1 .*      field2[:,1])
        @test (field2 .* field1)[:,1]   == (field1 .*      field2[:,1])
        @test (field2 .* field3)[:,1,1] == (field2[:,1] .* field3[:,1,1])
        @test (field1 .* field2_2)[:,1] == (field1 .*      field2_2[:,1])
        @test (field1 .* field2_2)[:,2] == (field1 .*      field2_2[:,2])
    end 
end

@testset "N-dimensional indexing" begin
    m, n, p = 2, 3, 4
    @testset for F in ( FullClenshawField,
                        FullGaussianField,
                        OctahedralGaussianField,
                        OctahedralClenshawField,
                        OctaminimalGaussianField,
                        HEALPixField,
                        OctaHEALPixField,
                        FullHEALPixField,
                        FullOctaHEALPixField,
                        )

        field = rand(F, m, n, p)
        @test field[:, 1, 1] isa F
        @test field[1] == field.data[1]
        @test field[1, 1, 1] == field.data[1, 1, 1]

        @test field[1:2, 1:2, 1:2] == field.data[1:2, 1:2, 1:2]
        @test field[1, 1, :] == field.data[1, 1, :]
        @test field[:, 1:2, 1:2].data == field.data[:, 1:2, 1:2]

        idx = CartesianIndex((1, 2, 3))
        @test field[idx] == field.data[idx]
        
        ids = CartesianIndices((m, n, p))
        @test field[ids] == field.data[ids]
        @test field[ids] isa Array
    end
end

@testset "Loop indexing" begin
    n = 2
    @testset for F in ( FullClenshawField,
                        FullGaussianField,
                        OctahedralGaussianField,
                        OctahedralClenshawField,
                        OctaminimalGaussianField,
                        HEALPixField,
                        OctaHEALPixField,
                        FullHEALPixField,
                        FullOctaHEALPixField,
                        )

        for s in ((n,), (n, n), (n, n, n), (n, n, n, n))
            field = zeros(F, s...)

            for k in eachlayer(field)
                for (j, ring) in enumerate(eachring(field))
                    for ij in ring
                        field[ij, k] = 1
                    end
                end
            end
            @test all(field .== 1)
        end
    end
end

@testset "AbstractField: GPU (JLArrays)" begin 
    NF = Float32
    @testset for F in ( 
        FullClenshawField,
        # FullGaussianField,            # don't test all for CI speedup
        OctahedralGaussianField,
        # OctahedralClenshawField,
        # OctaminimalGaussianField,
        HEALPixField,
        # OctaHEALPixField,
        # FullHEALPixField,
        # FullOctaHEALPixField,
    )
        s = (2, 3, 4)
        ndims = length(s)

        field_cpu = randn(F{NF}, s...)

        # constructors/adapt
        field = adapt(JLArray, field_cpu)
        field2 = Field(adapt(JLArray, field_cpu.data), field.grid)
        @test field == field2

        # broadcasting doesn't escape
        @test field  + field isa Field{NF, ndims, JLArray{NF, ndims}}
        @test field .+ field isa Field{NF, ndims, JLArray{NF, ndims}}
        @test field_cpu  + field_cpu isa Field{NF, ndims, Array{NF, ndims}}
        @test field_cpu .+ field_cpu isa Field{NF, ndims, Array{NF, ndims}}

        # getindex 
        @test field[1, :, :] isa JLArray{NF, 2}
        @test field[:, 1, 1] isa Field{NF, 1, JLArray{NF, 1}}
        for ring in eachring(field) 
            @test field[ring, :, :] == adapt(JLArray, field_cpu[ring, :, :])
            @test Array(field[ring, :, :]) == field_cpu[ring, :, :]
        end

        # setindex! 
        v = JLArray(rand(NF, s[3]))
        field[1, 1, :] .= v                # with .
        @test field[1, 1, :] == v

        v = JLArray(rand(NF, s[3]))
        field[2, 1, :] = v                 # without .
        @test field[2, 1, :] == v

        # with other field {Array}
        v = rand(F, s[1])
        field[:, 1, 1] = v                 # conversion to float64 -> float32
        @test Array(field[:, 1, 1].data) ≈ v

        # with other field {JLArray}
        v = adapt(JLArray, rand(F, s[1]))
        field[:, 1, 2] = v
        @test field[:, 1, 2].data ≈ v.data

        # fill 
        fill!(field, 2)
        @test all(field .== 2)
    end
end

@testset "Zonal mean" begin
    @testset for NF in (Int32, Int64, Float16, Float32, Float64)
        @testset for Grid in ( 
            FullClenshawGrid,
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
            FullHEALPixGrid,
            FullOctaHEALPixGrid,
        )
            nlat_half = 4
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)
            field = Field(1:npoints, grid)
        
            zm = zonal_mean(field)

            for (j, m) in enumerate(zonal_mean(field))
                @test m == sum(RingGrids.eachring(field)[j]) / RingGrids.get_nlon_per_ring(grid, j)
            end
        end
    end
end

@testset "nonparametric types" begin
    for M in (RingGrids, LowerTriangularArrays)
        @test M.nonparametric_type(Array) == Array
        @test M.nonparametric_type(Array{Float32}) == Array
        @test M.nonparametric_type(Array{Float32, 1}) == Array
        @test M.nonparametric_type(SubArray) == SubArray
        @test M.nonparametric_type(SubArray{Float32}) == SubArray
        @test M.nonparametric_type(SubArray{Float32, 1}) == SubArray
        @test M.nonparametric_type(SubArray{Float32, 1, Array}) == Array
        @test M.nonparametric_type(SubArray{Float32, 1, Array{Float32, 1}}) == Array
    end
end