@testset "Spherical distance on Earth" begin
    earth_circumference = 2π*RingGrids.DEFAULT_RADIUS
    for (point1, point2, expected) in [

        # longitude 
        ((0, 0), (1, 0), earth_circumference / 360),
        ((0, 0), (359, 0), earth_circumference / 360),
        ((0, 0), (180, 0), earth_circumference / 2),
        ((0, 0), (90, 0), earth_circumference / 4),

        # latitude
        ((0, 0), (0, 1), earth_circumference / 360),
        ((0, 0), (0, 10), 10*earth_circumference / 360),
        ((0, 0), (0, 90), earth_circumference / 4),
        ((0, 0), (0, -90), earth_circumference / 4),
    ]
        @test spherical_distance(point1, point2) ≈ expected
    end
end

@testset "Spherical distance invariants" begin

    # Ensure that function preserves the invariant that shifting in the longitude
    # direction does not change the distance.
    rand_lonlat() = (360 * rand(), 180 * rand() - 90)
    shift((lon, lat), lon_shift) = ((lon + lon_shift) % 360, lat)
    for _ in 1:100
        p1, p2 = rand_lonlat(), rand_lonlat()
        lon_shift = 360 * rand() 
        
        # Invariance in longitude shift.
        @test spherical_distance(p1, p2) ≈ spherical_distance(shift(p1, lon_shift), shift(p2, lon_shift))
        # Commutative.
        @test spherical_distance(p1, p2) ≈ spherical_distance(p2, p1)
        # Symmetry.
        @test spherical_distance(p1, p2) == spherical_distance((-p1[1], -p1[2]), (-p2[1], -p2[2]))
        # Identity.
        @test spherical_distance(p1, p1) == 0

        # Ensure that different ways to call the function are identical.
        @test spherical_distance(p1, p2) == spherical_distance(p1[1], p1[2], p2[1], p2[2])
        @test spherical_distance(p1, p2) == spherical_distance(p1, p2)
        @test spherical_distance(p1, p2) == spherical_distance(RingGrids.Haversine, p1, p2)
    end
end

@testset "Spherical distance with custom radius" begin
    custom_radius = 100
    circumference = 2 * π * custom_radius
    for (point1, point2, expected) in [
        ((0, 0), (1, 0), circumference / 360),
        ((0, 0), (359, 0), circumference / 360),
        ((0, 0), (180, 0), circumference / 2),
        ((0, 0), (90, 0), circumference / 4),

        ((0, 0), (0, 1), circumference / 360),
        ((0, 0), (0, 10), 10*circumference / 360),
        ((0, 0), (0, 90), circumference / 4),
        ((0, 0), (0, -90), circumference / 4),
    ]
        # There's a small numerical error, hence only approximate equality.
        @test spherical_distance(point1, point2, radius=custom_radius) ≈ expected
    end
end

@testset "Spherical distance type stability" begin
    for T in (Float16, Float32, Float64)
        p1, p2 = (rand(T), rand(T)), (rand(T), rand(T))
        @test typeof(spherical_distance(p1, p2)) == T

        # but if radius provided promote to that type
        for Tradius in (Float16, Float32, Float64)
            @test typeof(spherical_distance(p1, p2, radius=Tradius(1))) == promote_type(T, Tradius)
        end
    end

    # but integers always go to Float64 due to the trigonometric functions
    for T in (Int16, Int32, Int64)
        p1, p2 = (rand(T), rand(T)), (rand(T), rand(T))
        @test typeof(spherical_distance(p1, p2)) == Float64
    end
end
