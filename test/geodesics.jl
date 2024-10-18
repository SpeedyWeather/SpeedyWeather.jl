@testset "Haversine Earth" begin
    EARTH_CIRCUMFERENCE = 40075016.69  # Earth circumference in meters
    for (point1, point2, expected) in [
        ((0, 0), (1, 0), EARTH_CIRCUMFERENCE / 360),
        ((0, 0), (359, 0), EARTH_CIRCUMFERENCE / 360),
        ((0, 0), (180, 0), EARTH_CIRCUMFERENCE / 2),
        ((0, 0), (90, 0), EARTH_CIRCUMFERENCE / 4),
    ]
        @test RingGrids.Haversine(point1, point2) ≈ expected rtol=0.05
    end

    # Ensure that function preserves the invariant that shifting in the longitude
    # direction does not change the distance.
    rand_lonlat() = (360 * rand(), 180 * rand() - 90)
    shift((lon, lat), lon_shift) = ((lon + lon_shift) % 360, lat)
    for _ in 1:100
        p1, p2 = rand_lonlat(), rand_lonlat()
        lon_shift = 360 * rand() 
        
        # Invariance in longitude shift.
        @test RingGrids.Haversine(p1, p2) ≈ RingGrids.Haversine(shift(p1, lon_shift), shift(p2, lon_shift))
        # Commutative.
        @test RingGrids.Haversine(p1, p2) ≈ RingGrids.Haversine(p2, p1)
        # Symmetry.
        @test RingGrids.Haversine(p1, p2) == RingGrids.Haversine((-p1[1], -p1[2]), (-p2[1], -p2[2]))
        # Identity.
        @test RingGrids.Haversine(p1, p1) == 0

        # Ensure that different ways to call the function are identical.
        @test RingGrids.Haversine(p1, p2) == RingGrids.Haversine(p1[1], p1[2], p2[1], p2[2])
        @test RingGrids.Haversine(p1, p2) == RingGrids.spherical_distance(p1, p2)
        @test RingGrids.Haversine(p1, p2) == RingGrids.spherical_distance(RingGrids.Haversine, p1, p2)
    end
end

@testset "Haversine Custom Radius" begin
    custom_radius = 100
    circumference = 2 * π * custom_radius
    for (point1, point2, expected) in [
        ((0, 0), (1, 0), circumference / 360),
        ((0, 0), (359, 0), circumference / 360),
        ((0, 0), (180, 0), circumference / 2),
        ((0, 0), (90, 0), circumference / 4),
    ]
        # There's a small numerical error, hence only approximate equality.
        @test RingGrids.Haversine(point1, point2, radius=custom_radius) ≈ expected
    end
end