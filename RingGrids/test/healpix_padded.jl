using JLArrays

@testset "$Padded" for (Padded, Unpadded, nlat_halfs) in (
        (HEALPixPaddedGrid, HEALPixGrid, (16, 18, 24, 32)),
        (OctaHEALPixPaddedGrid, OctaHEALPixGrid, (8, 9, 16, 17)),
    )
    for nlat_half in nlat_halfs
        grid = Padded(nlat_half)
        healpix = Unpadded(nlat_half)
        nlat = RingGrids.get_nlat(grid)
        nlons = RingGrids.get_nlons(grid)
        nlons_healpix = RingGrids.get_nlons(healpix)

        @test all(nlons .>= 20)
        @test all(nlons .% nlons_healpix .== 0)
        @test [RingGrids.each_index_in_ring(Padded, j, nlat_half) for j in 1:nlat] == eachring(grid)     # nesting requires integer refinement
        @test all(iseven, nlons)
        @test nlons == reverse(nlons)                # north-south symmetric
        @test sum(nlons) == RingGrids.get_npoints(grid)
        @test RingGrids.get_nlat_half(Padded, RingGrids.get_npoints(grid)) == nlat_half
        @test RingGrids.get_latd(grid) == RingGrids.get_latd(healpix)
        @test sum(RingGrids.get_solid_angles(Padded, nlat_half) .* nlons) ≈ 4π

        # every HEALPix point is a padded point at k(i-1) + first
        for j in 1:nlat
            lond = RingGrids.get_lond_per_ring(Padded, nlat_half, j)
            lond_healpix = RingGrids.get_lond_per_ring(Unpadded, nlat_half, j)
            k, first = RingGrids.healpix_nesting(Padded, nlat_half, j)
            @test lond[k .* (0:(length(lond_healpix) - 1)) .+ first] ≈ lond_healpix
            @test lond[1] ≈ (RingGrids.hasoffset(Padded, nlat_half, j) ? 180 / nlons[j] : 0)
        end
    end

    # ring lengths from the pole as requested, identical to HEALPix from ring 5 on
    @test RingGrids.get_nlons(Padded(32))[1:8] == [20, 24, 24, 32, 20, 24, 28, 32]

    @testset "copyto! subsample and interpolate" begin
        for NF in (Float32, Float64), nlat_half in nlat_halfs[[1, end]], nlayers in ((), (3,))
            healpix = rand(NF, Unpadded(nlat_half), nlayers...)
            padded = zeros(NF, Padded(nlat_half), nlayers...)

            copyto!(padded, healpix)                    # interpolate
            healpix2 = zero(healpix)
            @test copyto!(healpix2, padded) === healpix2   # subsample
            @test healpix2 == healpix                   # nested points are copied exactly

            # interpolated values are bounded by the HEALPix values on that ring
            for (ring_healpix, ring_padded) in zip(eachring(healpix.grid), eachring(padded.grid))
                @test maximum(padded[ring_padded, :]) <= maximum(healpix[ring_healpix, :])
                @test minimum(padded[ring_padded, :]) >= minimum(healpix[ring_healpix, :])
            end

            # constant fields stay constant
            healpix .= 2
            copyto!(padded, healpix)
            @test all(padded .== 2)
        end

        # linear in longitude (away from the 0˚ wrap) is interpolated exactly
        nlat_half = nlat_halfs[1]
        healpix = Field(Float64, Unpadded(nlat_half))
        padded = Field(Float64, Padded(nlat_half))
        londs, _ = RingGrids.get_londlatds(healpix.grid)
        londs_padded, _ = RingGrids.get_londlatds(padded.grid)
        healpix.data .= londs
        copyto!(padded, healpix)
        for (j, ring) in enumerate(eachring(padded.grid))
            lond_healpix = RingGrids.get_lond_per_ring(Unpadded, nlat_half, j)
            inner = filter(ij -> lond_healpix[1] <= londs_padded[ij] <= lond_healpix[end], ring)
            @test padded[inner] ≈ londs_padded[inner]
        end

        @test_throws AssertionError copyto!(padded, zeros(Unpadded(nlat_halfs[end])))
        @test_throws DomainError Padded(nlat_halfs[1] - 2)

        # no allocations beyond the kernel launch
        healpix = rand(Unpadded(nlat_halfs[1]))
        padded = zeros(Padded(nlat_halfs[1]))
        copyto!(padded, healpix); copyto!(healpix, padded)
        @test @allocated(copyto!(padded, healpix)) < 2000
        @test @allocated(copyto!(healpix, padded)) < 2000

        # on a GPU-like architecture
        arch = architecture(JLArray)
        healpix_jl = on_architecture(arch, healpix)
        padded_jl = on_architecture(arch, zero(padded))
        copyto!(padded_jl, healpix_jl)
        @test Array(padded_jl.data) == padded.data
        healpix_jl .= 0
        copyto!(healpix_jl, padded_jl)
        @test Array(healpix_jl.data) == healpix.data
    end
end
