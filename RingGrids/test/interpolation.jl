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

@testset "Interpolate constant field, mixed precision and layers" begin
    # Float64 data with a (default) Float32 interpolator whose pole-average buffers fit all
    # layers: the pole averages must not be rounded to Float32, or constants are lost
    grid_in, grid_out = OctahedralGaussianGrid(8), FullGaussianGrid(12)
    nlayers = 3
    c = randn(Float64)
    A = fill!(Field(Float64, grid_in, nlayers), c)
    Aout = zeros(Float64, grid_out, nlayers)

    interpolator = RingGrids.interpolator(grid_out, grid_in; nlayers)
    @test eltype(interpolator.locator.north_pole_average) == Float32
    RingGrids.interpolate!(Aout, A, interpolator)
    @test all(a -> a ≈ c, Aout)
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

# `unsafe = true` only skips the checks on the input latitudes (e.g. for Reactant, which can't
# trace them, or for coordinates that are in range by construction), it never indexes out of bounds
@testset "update_locator! with and without safety checks (unsafe)" begin
    locator_arrays(L) = (L.js, L.ij_as, L.ij_bs, L.ij_cs, L.ij_ds, L.Δys, L.Δabs, L.Δcds)

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
            grid = Grid(4)
            (; nlat, npoints) = RingGrids.GridGeometry(grid; NF)

            # valid latitudes incl the poles, just inside the poles, on and between rings
            latd = NF.(RingGrids.get_latd(grid))
            θs = vcat(
                NF[90, prevfloat(NF(90)), 0, nextfloat(NF(-90)), -90],
                latd, (latd[1:(end - 1)] .+ latd[2:end]) / 2, 180 * rand(NF, 50) .- 90,
            )
            # any longitude is valid as it's periodic, incl. negative and ≥ 360˚
            λs = vcat(NF[0, 360, prevfloat(NF(360)), -eps(NF), 720], 720 * rand(NF, length(θs) - 5) .- 360)

            I_safe = RingGrids.interpolator(grid, length(θs); NF)
            I_unsafe = RingGrids.interpolator(grid, length(θs); NF)
            RingGrids.update_locator!(I_safe, λs, θs, unsafe = false)
            RingGrids.update_locator!(I_unsafe, λs, θs, unsafe = true)

            # identical for valid latitudes, i.e. the checks don't change the result
            @test all(map(==, locator_arrays(I_safe.locator), locator_arrays(I_unsafe.locator)))
            @test all(0 .<= I_safe.locator.js .<= nlat)
            @test all(0 .<= I_safe.locator.Δys .<= 1)

            # invalid latitudes: thrown with the safety checks
            θs_invalid = NF[90.5, -90.5, 1000, -Inf, NaN]
            λs_invalid = NF[10, 20, 30, 40, NaN]
            I = RingGrids.interpolator(grid, 1; NF)
            for (λ, θ) in zip(λs_invalid, θs_invalid)
                @test_throws ArgumentError RingGrids.update_locator!(I, [λ], [θ], unsafe = false)
            end

            # without the checks: north of 90˚N clamped to the north pole, south of -90˚N to
            # the south pole (extrapolated with Δy outside of [0, 1]), NaN to the north pole
            I = RingGrids.interpolator(grid, length(θs_invalid); NF)
            RingGrids.update_locator!(I, λs_invalid, θs_invalid, unsafe = true)
            (; js, ij_as, ij_bs, ij_cs, ij_ds, Δys) = I.locator
            @test js == [0, nlat, 0, nlat, 0]
            @test Δys[1] < 0 && Δys[2] > 1 && isnan(Δys[5])

            # but all indices stay in bounds (0 flags the north pole, -1 the south pole)
            @test all(ij -> ij == 0 || 1 <= ij <= npoints, vcat(ij_as, ij_bs))
            @test all(ij -> ij == -1 || 1 <= ij <= npoints, vcat(ij_cs, ij_ds))
            Aout = zeros(NF, length(θs_invalid))
            interpolate!(Aout, rand(NF, grid), I)       # no BoundsError (tests run with --check-bounds=yes)
            @test all(isfinite, Aout[1:3])
            @test isnan(Aout[5])
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

@testset "interpolate_2D!: batched matches per-layer" begin
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
            RingGrids.interpolate_2D!(field_out, field_in, interpolator)

            for k in 1:nlayers
                layer_in = zeros(NF, grid_in)
                layer_out = zeros(NF, grid_out)
                layer_in .= RingGrids.field_view(field_in, :, k)
                RingGrids.interpolate_2D!(layer_out, layer_in, interpolator)
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

@testset "interpolate_2D!: batched 4D field" begin
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
            RingGrids.interpolate_2D!(field_out, field_in, interpolator)

            for k in RingGrids.eachlayer(field_out, field_in)
                layer_in = zeros(NF, grid_in)
                layer_out = zeros(NF, grid_out)
                layer_in .= view(field_in.data, :, k)
                RingGrids.interpolate_2D!(layer_out, layer_in, interpolator)
                @test Array(layer_out) ≈ Array(view(field_out.data, :, k))
            end
        end
    end
end

@testset "interpolate_2D!: falls back for strided data" begin
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
    RingGrids.interpolate_2D!(out_strided, field_in, interpolator)
    RingGrids.interpolate_2D!(out_dense, dense_in, interpolator)

    @test Array(out_strided.data) ≈ Array(out_dense.data)
end

@testset "interpolate_2D!: constant field" begin
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
            RingGrids.interpolate_2D!(field_out, field_in, RingGrids.interpolator(grid_out, grid_in, NF = NF))

            for k in 1:nlayers
                @test all(Array(RingGrids.field_view(field_out, :, k)) .≈ NF(k))
            end
        end
    end
end

@testset "interpolate_2D!: contiguous views take the batched path" begin
    # a field wrapping a contiguous view (e.g. one time slice of a 4D array) is reshapeable
    # in O(1) and so is batched, not looped; the result must still match the dense one
    grid_in, grid_out = OctahedralGaussianGrid(8), OctahedralGaussianGrid(12)
    nlayers, ntime = 3, 2

    backing = randn(Float64, RingGrids.get_npoints(grid_in), nlayers, ntime)
    contiguous = view(backing, :, :, 2)
    @test RingGrids.is_flattenable(contiguous)
    @test !RingGrids.is_flattenable(view(backing, :, 1:2:nlayers, 1))

    field_in = Field(contiguous, grid_in)
    dense_in = Field(Array(contiguous), grid_in)

    out_view = zeros(Float64, grid_out, nlayers)
    out_dense = zeros(Float64, grid_out, nlayers)
    interpolator = RingGrids.interpolator(grid_out, grid_in, NF = Float64)
    RingGrids.interpolate_2D!(out_view, field_in, interpolator)
    RingGrids.interpolate_2D!(out_dense, dense_in, interpolator)

    @test Array(out_view.data) == Array(out_dense.data)
end

@testset "interpolate! forwards to interpolate_2D!" begin
    # without vertical positions interpolate! is a horizontal interpolation for every arity
    grid_in, grid_out = HEALPixGrid(8), FullGaussianGrid(12)
    interpolator = RingGrids.interpolator(grid_out, grid_in, NF = Float64)

    @testset for nlayers in (nothing, 4)
        field_in = isnothing(nlayers) ? randn(Float64, grid_in) : randn(Float64, grid_in, nlayers)
        out1 = isnothing(nlayers) ? zeros(Float64, grid_out) : zeros(Float64, grid_out, nlayers)
        out2 = deepcopy(out1)
        interpolate!(out1, field_in, interpolator)
        interpolate_2D!(out2, field_in, interpolator)
        @test out1.data == out2.data
        interpolate!(out1, field_in, interpolator.locator, interpolator.geometry)
        @test out1.data == out2.data
    end

    # vector output from a 2D field
    field_in = randn(Float64, grid_in)
    londs, latds = collect(-10.0:2.0:8.0), collect(0.0:5.0:45.0)
    I = RingGrids.interpolator(grid_in, length(londs), NF = Float64)
    RingGrids.update_locator!(I, londs, latds)
    v1, v2 = zeros(length(londs)), zeros(length(londs))
    interpolate!(v1, field_in, I)
    interpolate_2D!(v2, field_in, I)
    @test v1 == v2
    interpolate!(v1, field_in, I.locator, I.geometry)
    @test v1 == v2
end

@testset "interpolate_2D!: pole averages into the locator's buffers" begin
    # with pole-average buffers for enough layers the batched path writes into them instead of
    # allocating, a locator with fewer buffers (default nlayers = 1) falls back to allocating
    nlayers = 5
    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            OctaminimalGaussianGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
        )
        grid_in, grid_out = Grid(8), Grid(12)
        rings = grid_in.rings

        # data in a different number format than the interpolator, (Float32, Float64) uses the
        # buffers, (Float64, Float32) doesn't to keep the precision of the data
        @testset for (NF_data, NF_interp) in ((Float32, Float32), (Float64, Float64), (Float32, Float64), (Float64, Float32))
            field_in = randn(NF_data, grid_in, nlayers)
            out_buffers = zeros(NF_data, grid_out, nlayers)
            out_fallback = zeros(NF_data, grid_out, nlayers)

            I_buffers = RingGrids.interpolator(grid_out, grid_in; NF = NF_interp, nlayers)
            I_fallback = RingGrids.interpolator(grid_out, grid_in; NF = NF_interp)
            @test length(I_buffers.locator.north_pole_average) == nlayers
            @test length(I_fallback.locator.north_pole_average) == 1

            RingGrids.interpolate_2D!(out_buffers, field_in, I_buffers)
            RingGrids.interpolate_2D!(out_fallback, field_in, I_fallback)

            # ≈ as the buffers sum the pole rings in a different order than `mean`, and in the
            # interpolator's number format, which is the precision the comparison can expect
            @test isapprox(out_buffers.data, out_fallback.data, rtol = sqrt(eps(NF_interp)))

            if promote_type(NF_data, NF_interp) == NF_interp
                # buffers hold the pole averages, so the non-allocating path was taken
                @test I_buffers.locator.north_pole_average ≈ vec(Statistics.mean(field_in.data[rings[1], :], dims = 1))
                @test I_buffers.locator.south_pole_average ≈ vec(Statistics.mean(field_in.data[rings[end], :], dims = 1))
            else
                # buffers less precise than the data are left untouched, the allocating
                # fallback averages in the data's number format instead
                @test all(iszero, I_buffers.locator.north_pole_average)
            end

            # too small buffers are left untouched, the allocating fallback was taken
            @test all(iszero, I_fallback.locator.north_pole_average)
        end
    end

    # integer data takes the (rounding) fallback even with large enough buffers
    grid_in, grid_out = FullGaussianGrid(8), FullGaussianGrid(12)
    field_in = zeros(Int, grid_in, nlayers)
    field_in.data .= rand(1:10, size(field_in.data))
    out_buffers = zeros(Float64, grid_out, nlayers)
    out_fallback = zeros(Float64, grid_out, nlayers)
    I_buffers = RingGrids.interpolator(grid_out, grid_in; nlayers)
    RingGrids.interpolate_2D!(out_buffers, field_in, I_buffers)
    RingGrids.interpolate_2D!(out_fallback, field_in, RingGrids.interpolator(grid_out, grid_in))
    @test out_buffers.data == out_fallback.data
    @test all(iszero, I_buffers.locator.north_pole_average)
end

@testset "Interpolate ignoring NaN" begin
    @testset for NF in (Float32, Float64)
        grid_in, grid_out = FullGaussianGrid(16), OctahedralGaussianGrid(24)
        lond_in, _ = RingGrids.get_londlatds(grid_in)
        lond_out, latd_out = RingGrids.get_londlatds(grid_out)
        no_data = 0 .<= lond_in .< 90               # NaN "continent" in the first layer only

        A = randn(NF, grid_in, 2)
        A.data[no_data, 1] .= NaN
        plain = zeros(NF, grid_out, 2)
        interpolator = RingGrids.interpolator(plain, A)
        RingGrids.interpolate!(plain, A, interpolator)
        B = zeros(NF, grid_out, 2)
        RingGrids.interpolate_ignoring_nan!(B, A, interpolator)

        # unchanged wherever the plain interpolation is valid, incl. the layer without NaNs
        valid = .!isnan.(plain.data)
        @test B.data[valid] ≈ plain.data[valid]

        # NaN only without any valid neighbour, i.e. inside the continent but not next to it
        # (away from the poles, which interpolate from the average over the first/last ring)
        away_from_poles = abs.(latd_out) .< 80
        inside = away_from_poles .& (20 .< lond_out .< 70)
        outside = away_from_poles .& .!(0 .<= lond_out .< 90)
        @test all(isnan, B.data[inside, 1])
        @test any(isnan, plain.data[outside, 1])    # plain interpolation picks up NaNs next to the continent
        @test !any(isnan, B.data[outside, 1])       # these are filled from the valid neighbours

        # weights are renormalised over the valid neighbours: a constant stays constant
        c = randn(NF)
        A = fill!(Field(NF, grid_in, 2), c)
        A.data[no_data, 1] .= NaN
        RingGrids.interpolate_ignoring_nan!(B, A, interpolator)
        @test all(b -> b ≈ c, filter(!isnan, B.data))
    end
end
