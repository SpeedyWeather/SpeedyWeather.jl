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

            # check unmasked_indices returns exactly the positions where mask is false
            @test issorted(indices)
            @test indices == sort(findall(.~mask.data))

            # 2D (single-layer) round-trip
            src2d = rand(NF, grid)
            src2d_copy = copy(src2d)
            array2d = zeros(NF, n_unmasked)
            copy_unmasked!(array2d, src2d, indices)

            @test all(array2d[i] == src2d[indices[i]] for i in 1:n_unmasked)

            # overwrite the unmasked positions in src, then copy back and check restoration
            src2d .= zero(NF)
            copy_unmasked!(src2d, array2d, indices)

            @test all(src2d[indices[i]] ≈ src2d_copy[indices[i]] for i in 1:n_unmasked)
            # masked positions were not touched (still zero)
            masked_ij = findall(mask.data)
            @test all(src2d[ij] == zero(NF) for ij in masked_ij)

            # 3D (multi-layer) round-trip
            nlayers = 5
            src3d = rand(NF, grid, nlayers)
            src3d_copy = copy(src3d)
            array3d = zeros(NF, n_unmasked, nlayers)
            copy_unmasked!(array3d, src3d, indices)

            @test all(array3d[i, k] == src3d[indices[i], k] for i in 1:n_unmasked, k in 1:nlayers)

            src3d .= zero(NF)
            copy_unmasked!(src3d, array3d, indices)

            @test all(src3d[indices[i], k] ≈ src3d_copy[indices[i], k] for i in 1:n_unmasked, k in 1:nlayers)
            @test all(src3d[ij, k] == zero(NF) for ij in masked_ij, k in 1:nlayers)
        end
    end
end
