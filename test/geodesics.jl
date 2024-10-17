@testset "Haversine" begin
    for (point1, point2, expected) in [
        # Check that distance at Equator is roughly 111321.0 m
        ((0, 0), (1, 0), 111321.0),
        ((0, 0), (359, 0), 111321.0),
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
        
        @test RingGrids.Haversine(p1, p2) ≈ RingGrids.Haversine(shift(p1, lon_shift), shift(p2, lon_shift))
        @test RingGrids.Haversine(p1, p1) == 0

        # Ensure that different ways to call the function are identical.
        @test RingGrids.Haversine(p1, p2) == RingGrids.Haversine(p1[1], p1[2], p2[1], p2[2])
        @test RingGrids.Haversine(p1, p2) == RingGrids.spherical_distance(p1, p2)
        @test RingGrids.Haversine(p1, p2) == RingGrids.spherical_distance(RingGrids.Haversine, p1, p2)
    end
end