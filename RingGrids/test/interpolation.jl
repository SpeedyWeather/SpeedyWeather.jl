import Random
import Statistics

@testset "Interpolate constant field" begin
    npoints = 100

    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )

        @testset for NF in (Float32, Float64)

            grid = Grid(8)
            A = randn(NF, grid)             # a very low resolution grid
            c = randn(NF)
            A .= c                          # constant data globally

            λs = 360 * rand(NF, npoints)          # some longitudes in [0˚, 360˚E]
            θs = 180 * rand(NF, npoints) .- 90    # some latitudes in [-90˚, 90˚N]
            As = RingGrids.interpolate(λs, θs, A)

            for a in As
                @test a ≈ c
            end
        end
    end
end

@testset "Interpolate zonally-constant field" begin
    npoints = 1000

    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )

        @testset for NF in (Float32, Float64)

            grid = Grid(32)
            A = zeros(NF, grid)                 # some resolution
            G = RingGrids.GridGeometry(A)
            lat1 = G.latd[2]                    # latitude of first ring

            for (j, ring) in enumerate(RingGrids.eachring(A))
                θ = G.latd[j + 1]     # G.latd also includes 90˚N hence +1
                for ij in ring
                    A[ij] = θ
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            λs = 360 * rand(NF, npoints)              # some longitudes in [0˚, 360˚E]
            θs = 2lat1 * rand(NF, npoints) .- lat1    # some latitudes in [-90˚, 90˚N]

            As = RingGrids.interpolate(λs, θs, A; NF)

            for (a, θ) in zip(As, θs)
                @test a ≈ θ
            end
        end
    end
end

@testset "Interpolate meridionally-constant field" begin
    npoints = 1000

    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )

        @testset for NF in (Float32, Float64)

            grid = Grid(32)
            A = zeros(NF, grid)             # some resolution
            G = RingGrids.GridGeometry(A)
            lat1 = G.latd[2]                # latitude of first ring

            # TEST FROM 60˚ to 300˚ to not interpolate across 0/360˚E
            # where this test won't work because lon have a sharp jump
            # and aren't linear across the prime meridian
            # but that differently further down
            for (j, ring) in enumerate(RingGrids.eachring(A))
                for ij in ring
                    A[ij] = G.londs[ij]
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            λs = 240 * rand(NF, npoints) .+ 60        # some longitudes in [60˚, 300˚E]
            θs = 2lat1 * rand(NF, npoints) .- lat1    # some latitudes in (-90˚, 90˚N)

            As = RingGrids.interpolate(λs, θs, A; NF)

            for (a, λ) in zip(As, λs)
                @test a ≈ λ rtol = 1.0e-3
            end

            f(λ) = λ > 180 ? λ - 360 : λ          # 0-360˚ to -180˚-180˚E

            # TEST FROM -120˚ to 120˚ to still test the indexing across the
            # prime meridian
            for (j, ring) in enumerate(RingGrids.eachring(A))
                for ij in ring
                    A[ij] = f(G.londs[ij])
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            λs = 240 * rand(NF, npoints) .- 120       # some longitudes in [-120˚, 120˚E]
            θs = 2lat1 * rand(NF, npoints) .- lat1    # some latitudes in (-90˚, 90˚N)

            As = RingGrids.interpolate(λs, θs, A; NF)

            for (a, λ) in zip(As, λs)
                @test a ≈ λ rtol = 2.5e-3
            end
        end
    end
end

@testset "Find latitude rings and weights" begin
    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )

        @testset for nlat_half in [4, 8, 16]
            grid = Grid(nlat_half)
            G = RingGrids.GridGeometry(grid)
            latd = G.latd

            n = length(latd) - 1
            Δs = rand(n)

            r = Random.randperm(n)
            θs = latd[1:(end - 1)] .+ diff(latd) .* Δs

            θs = θs[r]
            Δs = Δs[r]
            js, Δys = RingGrids.find_rings(θs, latd)

            for (i, (j, Δref, Δ)) in enumerate(zip(js, Δs, Δys))
                @test j == r[i] - 1
                @test Δref ≈ Δ
            end
        end
    end
end

@testset "Interpolate between grids" begin
    @testset for NF in (Float32, Float64)
        @testset for Grid in (
                FullGaussianGrid,
                FullClenshawGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                OctaminimalGaussianGrid,
                HEALPixGrid,
                OctaHEALPixGrid,
            )

            # create some gridded field
            nlat_half_src = 16
            A = randn(NF, Grid(nlat_half_src))

            # interpolate to FullGaussianGrid and back and compare
            nlat_half = 32
            grid = FullGaussianGrid(nlat_half)
            A_interpolated = RingGrids.interpolate(grid, A)
            A2 = zero(A)
            RingGrids.interpolate!(A2, A_interpolated)

            # just check that it's not completely off
            @test A ≈ A2 rtol = 5.0e-1 atol = 5.0e-1
        end
    end
end

@testset "3/4D interpolation interfaces" begin
    A = randn(OctahedralGaussianField, 16, 2)
    B = zeros(FullGaussianField, 16, 2)
    C = zeros(FullGaussianField, 16, 2)

    RingGrids.interpolate!(B, A)

    interpolator = RingGrids.interpolator(C, A)
    RingGrids.interpolate!(C, A, interpolator)

    @test B == C

    A = randn(OctahedralGaussianField, 8, 3, 2)
    B = zeros(FullGaussianField, 8, 3, 2)
    C = zeros(FullGaussianField, 8, 3, 2)

    RingGrids.interpolate!(B, A)

    interpolator = RingGrids.interpolator(C, A)
    RingGrids.interpolate!(C, A, interpolator)

    @test B == C
end

@testset "Grid cell average" begin
    for Grid in (
            FullGaussianGrid,
            FullClenshawGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )
        for NF in (Float32, Float64)

            for nlat_half in (8, 16)
                full_grid = RingGrids.full_grid_type(Grid)(nlat_half)
                full_field = randn(NF, full_grid) .+ 3

                for nlat_half in (4, 8)
                    grid = Grid(nlat_half)
                    field = zeros(Float32, grid)
                    RingGrids.grid_cell_average!(field, full_field)
                    @test Statistics.std(field) < Statistics.std(full_field)
                    @test minimum(field) >= minimum(full_field)
                    @test maximum(field) <= maximum(full_field)
                end
            end
        end
    end
end

@testset "Batched interpolation matches per-layer" begin
    # A 3D field interpolates all its layers in one launch; the result must agree exactly
    # with interpolating each layer separately as a 2D field.
    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )
        @testset for NF in (Float32, Float64)
            grid_in = Grid(8)
            grid_out = Grid(12)
            nlayers = 5

            field_in = randn(NF, grid_in, nlayers)
            field_out = zeros(NF, grid_out, nlayers)

            interpolator = RingGrids.interpolator(grid_out, grid_in, NF = NF)
            RingGrids.interpolate!(field_out, field_in, interpolator)

            for k in 1:nlayers
                layer_in = zeros(NF, grid_in)
                layer_out = zeros(NF, grid_out)
                layer_in .= RingGrids.field_view(field_in, :, k)
                RingGrids.interpolate!(layer_out, layer_in, interpolator)
                # ≈ rather than ==: the batched path averages the pole rings with a
                # reduction over all layers at once, which accumulates in a different
                # order than the single-layer reduction. On some grids that differs in
                # the last bit, and the difference reaches only those output points that
                # take a pole value (observed: HEALPix-family grids in Float64, 1 ULP).
                @test Array(layer_out) ≈ Array(RingGrids.field_view(field_out, :, k))
            end
        end
    end
end

@testset "Batched interpolation of a 4D field" begin
    # trailing dimensions are collapsed into one for the batched launch; the result must
    # match interpolating each (layer, trailing) slice on its own
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        @testset for NF in (Float32, Float64)
            grid_in = Grid(8)
            grid_out = Grid(12)
            n2, n3 = 4, 3

            field_in = randn(NF, grid_in, n2, n3)
            field_out = zeros(NF, grid_out, n2, n3)

            interpolator = RingGrids.interpolator(grid_out, grid_in, NF = NF)
            RingGrids.interpolate!(field_out, field_in, interpolator)

            for k in RingGrids.eachlayer(field_out, field_in)
                layer_in = zeros(NF, grid_in)
                layer_out = zeros(NF, grid_out)
                layer_in .= view(field_in.data, :, k)
                RingGrids.interpolate!(layer_out, layer_in, interpolator)
                @test Array(layer_out) ≈ Array(view(field_out.data, :, k))
            end
        end
    end
end

@testset "Batched interpolation falls back for non-dense data" begin
    # a field whose data is a non-contiguous view cannot be reshaped in O(1), so it takes
    # the per-layer loop instead; either way the answer must be the same
    grid_in, grid_out = OctahedralGaussianGrid(8), OctahedralGaussianGrid(12)
    nlayers = 3

    backing = randn(Float64, RingGrids.get_npoints(grid_in), 2 * nlayers)
    strided = view(backing, :, 1:2:(2 * nlayers))        # every other column
    @test !(strided isa DenseArray)

    field_in = Field(strided, grid_in)
    dense_in = Field(Array(strided), grid_in)

    out_strided = zeros(Float64, grid_out, nlayers)
    out_dense = zeros(Float64, grid_out, nlayers)
    interpolator = RingGrids.interpolator(grid_out, grid_in, NF = Float64)
    RingGrids.interpolate!(out_strided, field_in, interpolator)
    RingGrids.interpolate!(out_dense, dense_in, interpolator)

    @test Array(out_strided.data) ≈ Array(out_dense.data)
end

@testset "Batched interpolation of a constant field" begin
    # constants must survive, including at the poles where the batched kernel uses
    # per-layer pole averages rather than the single-layer kernel's scalars
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        @testset for NF in (Float32, Float64)
            grid_in = Grid(8)
            grid_out = Grid(16)
            nlayers = 3

            field_in = zeros(NF, grid_in, nlayers)
            for k in 1:nlayers
                RingGrids.field_view(field_in, :, k) .= NF(k)
            end

            field_out = zeros(NF, grid_out, nlayers)
            RingGrids.interpolate!(field_out, field_in, RingGrids.interpolator(grid_out, grid_in, NF = NF))

            for k in 1:nlayers
                @test all(Array(RingGrids.field_view(field_out, :, k)) .≈ NF(k))
            end
        end
    end
end
