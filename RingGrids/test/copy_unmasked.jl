@testset "copy_unmasked! round-trip" begin
    for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            HEALPixGrid,
        )
        for NF in (Float32, Float64)
            grid = Grid(8)
            npoints = RingGrids.get_npoints(grid)

            # random boolean mask: ~half the points masked (true = masked, false = unmasked)
            mask = rand(Bool, grid)
            indices = unmasked_indices(mask)
            n_unmasked = length(indices)
            @test n_unmasked == length(mask) - sum(mask)

            # check unmasked_indices returns exactly the positions where mask is false
            @test issorted(indices)
            @test indices == sort(findall(.~mask))

            # 2D (single-layer) round-trip
            field2d = rand(NF, grid)
            field2d_copy = copy(field2d)
            array2d = zeros(NF, n_unmasked)
            copy_unmasked!(array2d, field2d, indices)
            @test array2d == field2d[.~mask]

            # overwrite the unmasked positions in src, then copy back and check restoration
            field2d .= zero(NF)
            copy_unmasked!(field2d, array2d, indices)
            
            # agreement in non-masked elements
            @test field2d[.~mask] == field2d_copy[.~mask]
            
            # masked positions were not touched (still zero)
            @test all(field2d[mask] .== 0)

            # 3D (multi-layer) round-trip
            nlayers = 5
            field3d = rand(NF, grid, nlayers)
            field3d_copy = copy(field3d)
            array3d = zeros(NF, n_unmasked, nlayers)
            copy_unmasked!(array3d, field3d, indices)
            @test array3d == field3d[.~mask, :]

            # copy back
            field3d .= zero(NF)
            copy_unmasked!(field3d, array3d, indices)

            # agreement in non-masked elements
            @test field3d[.~mask, :] == field3d_copy[.~mask, :]
            
            # masked positions were not touched (still zero)
            @test all(field3d[mask, :] .== 0)
        end
    end
end
