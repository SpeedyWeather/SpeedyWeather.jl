@testset "HEALPixPaddedGrid" begin

    @testset "ring structure follows HEALPix with padded caps" begin
        for nlat_half in (8, 16, 32, 144)
            pad = RingGrids.npoints_padding(HEALPixPaddedGrid)
            nlat = RingGrids.get_nlat(HEALPixPaddedGrid, nlat_half)
            @test nlat == RingGrids.get_nlat(HEALPixGrid, nlat_half)

            for j in 1:nlat
                nlon_padded = RingGrids.get_nlon_per_ring(HEALPixPaddedGrid, nlat_half, j)
                nlon_healpix = RingGrids.get_nlon_per_ring(HEALPixGrid, nlat_half, j)

                # never fewer points than HEALPix, never more than the belt
                @test nlon_padded >= nlon_healpix
                @test nlon_padded <= 2nlat_half
                @test nlon_padded == min(nlon_healpix + pad, 2nlat_half)

                # a real FFT needs an even ring length
                @test iseven(nlon_padded)
            end

            # the belt itself is untouched
            @test RingGrids.get_nlon_max(HEALPixPaddedGrid, nlat_half) == 2nlat_half
        end
    end

    @testset "latitudes are HEALPix's, unchanged" begin
        for nlat_half in (8, 32, 144)
            @test RingGrids.get_latd(HEALPixPaddedGrid, nlat_half) ==
                RingGrids.get_latd(HEALPixGrid, nlat_half)
        end
    end

    @testset "npoints and its inverse are consistent" begin
        for nlat_half in (2, 8, 16, 32, 144)
            npoints = RingGrids.get_npoints(HEALPixPaddedGrid, nlat_half)
            @test npoints == sum(
                RingGrids.get_nlon_per_ring(HEALPixPaddedGrid, nlat_half, j)
                    for j in 1:RingGrids.get_nlat(HEALPixPaddedGrid, nlat_half)
            )
            @test npoints >= RingGrids.get_npoints(HEALPixGrid, nlat_half)
            @test RingGrids.get_nlat_half(HEALPixPaddedGrid, npoints) == nlat_half
        end
    end

    @testset "ring indices partition the grid" begin
        for nlat_half in (8, 16, 32)
            grid = HEALPixPaddedGrid(nlat_half)
            rings = RingGrids.eachring(grid)
            nlat = RingGrids.get_nlat(grid)
            @test length(rings) == nlat
            @test first(first(rings)) == 1
            @test last(last(rings)) == RingGrids.get_npoints(HEALPixPaddedGrid, nlat_half)
            for j in 1:nlat
                @test length(rings[j]) == RingGrids.get_nlon_per_ring(grid, j)
                # each_index_in_ring must agree with the precomputed ranges
                @test rings[j] == RingGrids.each_index_in_ring(HEALPixPaddedGrid, j, nlat_half)
                j > 1 && @test first(rings[j]) == last(rings[j - 1]) + 1
            end
        end
    end

    @testset "quadrature weights integrate the sphere" begin
        # The ring *areas* are HEALPix's — same latitudes, so the same colatitude bands — but they
        # are shared among more points on the cap rings, so the per-point solid angle differs.
        for nlat_half in (8, 32, 144)
            weights = RingGrids.get_quadrature_weights(HEALPixPaddedGrid, nlat_half)
            # HEALPixGrid defines get_solid_angles directly and has no get_quadrature_weights,
            # so compare against the equal-area ring weights it is built from
            @test weights == RingGrids.equal_area_weights(HEALPixGrid, nlat_half)
            @test sum(weights) ≈ 2                       # ∫₀^π sin θ dθ = 2

            solid_angles = RingGrids.get_solid_angles(HEALPixPaddedGrid, nlat_half)
            nlons = [
                RingGrids.get_nlon_per_ring(HEALPixPaddedGrid, nlat_half, j)
                    for j in 1:RingGrids.get_nlat(HEALPixPaddedGrid, nlat_half)
            ]
            @test sum(nlons .* solid_angles) ≈ 4π        # the whole sphere

            # Unlike HEALPixGrid this is NOT equal area: a padded cap ring shares HEALPix's ring
            # area among more points, so its pixels are smaller than the belt's. Only meaningful
            # once the grid is fine enough that some cap ring is not saturated to the belt width.
            if nlat_half > 8
                @test !all(≈(first(solid_angles)), solid_angles)
                @test minimum(solid_angles) < maximum(solid_angles)
            end
        end
    end

    @testset "fields work and integrate correctly" begin
        for nlat_half in (8, 16)
            grid = HEALPixPaddedGrid(nlat_half)
            field = ones(Float64, grid)
            @test length(field) == RingGrids.get_npoints(HEALPixPaddedGrid, nlat_half)

            # the area-weighted mean of a constant field is that constant
            solid_angles = RingGrids.get_solid_angles(HEALPixPaddedGrid, nlat_half)
            total = sum(
                sum(field[r]) * solid_angles[j]
                    for (j, r) in enumerate(RingGrids.eachring(grid))
            )
            @test total ≈ 4π
        end
    end

    @testset "longitude offsets follow the HEALPix convention" begin
        for nlat_half in (8, 32)
            nlat = RingGrids.get_nlat(HEALPixPaddedGrid, nlat_half)
            for j in 1:nlat
                lond = RingGrids.get_lond_per_ring(HEALPixPaddedGrid, nlat_half, j)
                nlon = RingGrids.get_nlon_per_ring(HEALPixPaddedGrid, nlat_half, j)
                @test length(lond) == nlon
                @test issorted(lond)
                @test all(0 .<= lond .< 360)
                # equidistant
                @test all(≈(360 / nlon), diff(lond))
                # same offset flag as HEALPixGrid on the same ring
                @test RingGrids.hasoffset(HEALPixPaddedGrid, nlat_half, j) ==
                    RingGrids.hasoffset(HEALPixGrid, nlat_half, j)
            end
            # asking for a whole-grid offset is ambiguous here, as for HEALPixGrid
            @test_throws ArgumentError RingGrids.hasoffset(HEALPixPaddedGrid)
        end
    end
end
