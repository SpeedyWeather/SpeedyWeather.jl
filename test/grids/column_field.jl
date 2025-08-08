using JLArrays
using Adapt

RINGGRIDS_DEFAULT_NF = SpeedyWeather.RingGrids.DEFAULT_NF

@testset "ColumnField types" begin
    # Test type hierarchy and properties
    @test ColumnField <: AbstractField
    @test FullColumnField <: ColumnField
    @test ReducedColumnField <: ColumnField
    
    # Test type aliases
    @test ColumnField2D == ColumnField{T, 1} where T
    @test ColumnField3D == ColumnField{T, 2} where T  
    @test ColumnField4D == ColumnField{T, 3} where T
end

@testset "ColumnField constructors" begin
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
            nlat_half = 4
            nlayers = 10
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)

            # Test basic constructor with data and grid
            data = rand(NF, nlayers, npoints)
            field = ColumnField(data, grid)
            @test field.data === data
            @test field.grid === grid
            @test eltype(field) == NF
            @test size(field) == (nlayers, npoints)

            # Test dimension mismatch error
            wrong_data = rand(NF, nlayers, npoints - 1)
            @test_throws DimensionMismatch ColumnField(wrong_data, grid)

            # Test default constructors
            field1 = ColumnField(grid, nlayers)
            field2 = ColumnField(NF, grid, nlayers)
            @test size(field1) == (nlayers, npoints)
            @test size(field2) == (nlayers, npoints)
            @test eltype(field1) == RINGGRIDS_DEFAULT_NF
            @test eltype(field2) == NF
        end
    end
end

@testset "ColumnField undef constructors" begin
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
            nlat_half = 4
            nlayers = 8
            npoints = RingGrids.get_npoints(Grid, nlat_half)

            # Get the corresponding ColumnField type
            F = RingGrids.field_type(Grid)
            ColumnF = ColumnField{NF, 2, Array{NF, 2}, Grid}

            # Test various undef constructor patterns
            field1 = ColumnF(undef, nlayers, nlat_half)
            @test size(field1) == (nlayers, npoints)
            @test eltype(field1) == NF
            @test field1.grid isa Grid

            # Test with additional dimensions
            field2 = ColumnF(undef, nlayers, nlat_half, 3)
            @test size(field2) == (nlayers, npoints, 3)
            @test eltype(field2) == NF

            # Test generic ColumnField constructor
            field3 = ColumnField(undef, nlayers, nlat_half)
            @test size(field3) == (nlayers, npoints)
            @test eltype(field3) == RINGGRIDS_DEFAULT_NF
        end
    end
end

@testset "ColumnField indexing" begin
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
            nlat_half = 4
            nlayers = 6
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)

            field = ColumnField(rand(NF, nlayers, npoints), grid)

            # Test getindex
            @test field[1, 1] isa NF
            @test field[1:2, 1] isa Vector{NF}
            @test field[1, 1:3] isa Vector{NF}
            @test size(field[1:3, 1:5]) == (3, 5)

            # Test setindex
            original_val = field[1, 1]
            field[1, 1] = NF(42)
            @test field[1, 1] == NF(42)
            field[1, 1] = original_val  # restore

            # Test linear indexing
            for i in eachindex(field)
                field[i] = NF(i)
            end
            @test field[1] == NF(1)
            @test field[end] == NF(length(field))
        end
    end
end

@testset "ColumnField transpose operations" begin
    for NF in (Float32, Float64)
        for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            nlat_half = 4
            nlayers = 5
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)

            # Create a regular Field
            field_data = rand(NF, npoints, nlayers)
            field = Field(field_data, grid)

            # Test transpose from Field to ColumnField
            column_field = transpose(field)
            @test column_field isa ColumnField
            @test size(column_field) == (nlayers, npoints)
            @test eltype(column_field) == NF

            # Test transpose back from ColumnField to Field
            field_back = transpose(column_field)
            @test field_back isa Field
            @test size(field_back) == (npoints, nlayers)
            @test eltype(field_back) == NF

            # Test unsafe transpose operations
            field_copy = Field(copy(field_data), grid)
            column_field_unsafe = transpose!(field_copy)
            @test column_field_unsafe isa ColumnField
            @test size(column_field_unsafe) == (nlayers, npoints)

            # Create ColumnField directly
            column_data = rand(NF, nlayers, npoints)
            column_field = ColumnField(column_data, grid)
            
            # Test unsafe transpose back
            column_copy = ColumnField(copy(column_data), grid)
            field_unsafe = transpose!(column_copy)
            @test field_unsafe isa Field
            @test size(field_unsafe) == (npoints, nlayers)
        end
    end
end

@testset "ColumnField similar operations" begin
    for NF in (Float32, Float64)
        for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            nlat_half = 4
            nlayers = 6
            grid = Grid(nlat_half)
            
            field = ColumnField(rand(NF, nlayers, RingGrids.get_npoints(grid)), grid)

            # Test similar with same type, new size
            new_nlat_half = 8
            new_nlayers = 10
            similar_field = similar(field, new_nlayers, new_nlat_half)
            @test similar_field isa ColumnField
            @test eltype(similar_field) == NF
            @test size(similar_field, 1) == new_nlayers
            @test size(similar_field, 2) == RingGrids.get_npoints(Grid, new_nlat_half)

            # Test similar with new type and size
            similar_field_f32 = similar(field, Float32, new_nlayers, new_nlat_half)
            @test similar_field_f32 isa ColumnField
            @test eltype(similar_field_f32) == Float32
            @test size(similar_field_f32, 1) == new_nlayers
            @test size(similar_field_f32, 2) == RingGrids.get_npoints(Grid, new_nlat_half)

            # Test similar with additional dimensions
            similar_field_3d = similar(field, new_nlayers, new_nlat_half, 3)
            @test size(similar_field_3d) == (new_nlayers, RingGrids.get_npoints(Grid, new_nlat_half), 3)
        end
    end
end

@testset "ColumnField broadcasting" begin
    for NF in (Float32, Float64)
        for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            nlat_half = 4
            nlayers = 5
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)

            field1 = ColumnField(rand(NF, nlayers, npoints), grid)
            field2 = ColumnField(rand(NF, nlayers, npoints), grid)

            # Test basic arithmetic broadcasting
            result_add = field1 .+ field2
            @test result_add isa ColumnField
            @test size(result_add) == size(field1)
            @test eltype(result_add) == NF

            result_mul = field1 .* field2
            @test result_mul isa ColumnField
            @test size(result_mul) == size(field1)

            # Test broadcasting with scalars
            result_scalar = field1 .+ NF(2)
            @test result_scalar isa ColumnField
            @test size(result_scalar) == size(field1)

            result_scalar2 = NF(3) .* field1
            @test result_scalar2 isa ColumnField
            @test size(result_scalar2) == size(field1)

            # Test in-place broadcasting
            field_copy = ColumnField(copy(field1.data), grid)
            field_copy .+= field2
            @test field_copy isa ColumnField
            @test size(field_copy) == size(field1)
        end
    end
end

@testset "ColumnField Array conversion" begin
    for Grid in (FullGaussianGrid, FullClenshawGrid)  # Only test full grids
        nlat_half = 4
        nlayers = 6
        grid = Grid(nlat_half)
        npoints = RingGrids.get_npoints(grid)
        nlat = RingGrids.get_nlat(grid)

        field = ColumnField(rand(Float64, nlayers, npoints), grid)

        # Test Array conversion for full grids
        if field isa FullColumnField
            matrix = Array(field, Matrix)
            @test matrix isa Matrix
            @test size(matrix, 1) == nlat
            @test size(matrix, 2) == nlayers
        end
    end
end

@testset "ColumnField arithmetic with Field" begin
    for NF in (Float32, Float64)
        for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            nlat_half = 4
            nlayers = 5
            grid = Grid(nlat_half)
            npoints = RingGrids.get_npoints(grid)

            # Create Field and ColumnField
            field = Field(rand(NF, npoints, nlayers), grid)
            column_field = ColumnField(rand(NF, nlayers, npoints), grid)

            # Test add! operations
            field_copy = Field(copy(field.data), grid)
            add!(field_copy, column_field)
            @test field_copy isa Field
            @test size(field_copy) == size(field)

            column_field_copy = ColumnField(copy(column_field.data), grid)
            add!(column_field_copy, field)
            @test column_field_copy isa ColumnField
            @test size(column_field_copy) == size(column_field)
        end
    end
end

@testset "ColumnField: GPU (JLArrays)" begin
    NF = Float32
    @testset for Grid in (
        FullClenshawGrid,
        OctahedralGaussianGrid,
        HEALPixGrid,
    )
        nlat_half = 4
        nlayers = 6
        grid = Grid(nlat_half)
        npoints = RingGrids.get_npoints(grid)

        # Create CPU ColumnField
        column_field_cpu = ColumnField(rand(NF, nlayers, npoints), grid)

        # Test adapt to GPU
        column_field_gpu = adapt(JLArray, column_field_cpu)
        @test column_field_gpu isa ColumnField{NF, 2, JLArray{NF, 2}}
        @test size(column_field_gpu) == size(column_field_cpu)

        # Test constructor with GPU array
        gpu_data = JLArray(rand(NF, nlayers, npoints))
        column_field_gpu2 = ColumnField(gpu_data, grid)
        @test column_field_gpu2 isa ColumnField{NF, 2, JLArray{NF, 2}}

        # Test broadcasting preserves GPU type
        result = column_field_gpu .+ column_field_gpu
        @test result isa ColumnField{NF, 2, JLArray{NF, 2}}

        # Test indexing returns GPU arrays
        @test column_field_gpu[1, :] isa JLArray{NF, 1}
        @test column_field_gpu[1:2, 1:3] isa JLArray{NF, 2}

        # Test setindex! with GPU arrays
        gpu_vec = JLArray(rand(NF, npoints))
        column_field_gpu[1, :] = gpu_vec
        @test column_field_gpu[1, :] ≈ gpu_vec

        # Test fill!
        fill!(column_field_gpu, NF(42))
        @test all(Array(column_field_gpu) .== NF(42))
    end
end

@testset "ColumnField type utilities" begin
    # Test nonparametric_type
    @test SpeedyWeather.nonparametric_type(ColumnField) == ColumnField
    @test SpeedyWeather.nonparametric_type(ColumnField{Float32}) == ColumnField
    @test SpeedyWeather.nonparametric_type(ColumnField{Float32, 2}) == ColumnField

    # Test grid_type extraction
    grid = FullGaussianGrid(4)
    field = ColumnField(rand(5, RingGrids.get_npoints(grid)), grid)
    @test RingGrids.grid_type(typeof(field)) == typeof(grid)

    # Test array_type extraction
    @test Architectures.array_type(typeof(field)) == Array{Float64, 2}
end

@testset "ColumnField error handling" begin
    grid = FullGaussianGrid(4)
    npoints = RingGrids.get_npoints(grid)

    # Test dimension mismatch in constructor
    wrong_data = rand(5, npoints + 1)  # Wrong number of points
    @test_throws DimensionMismatch ColumnField(wrong_data, grid)

    # Test bounds errors in transpose operations
    field = ColumnField(rand(5, npoints), grid)
    wrong_scratch = rand(3, 4)  # Wrong size scratch array
    @test_throws BoundsError RingGrids.transpose_unsafe!(field, wrong_scratch)
end
