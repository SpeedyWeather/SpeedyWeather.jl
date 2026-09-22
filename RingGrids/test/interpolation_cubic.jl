@testset "CubicInterpolator: exact on the grid points" begin
    for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        grid = Grid(24)
        npoints = RingGrids.get_npoints(grid)
        londs, latds = RingGrids.get_londlatds(grid)

        I = CubicInterpolator(grid, npoints; NF = Float64)
        RingGrids.update_locator!(I, londs, latds)

        # interpolating onto the grid points themselves must return the field unchanged
        for f in ((lo, la) -> 1.0, (lo, la) -> sind(la), (lo, la) -> cosd(la)^4 * cosd(4lo))
            A = Field([f(londs[i], latds[i]) for i in 1:npoints], grid)
            out = zeros(Float64, npoints)
            RingGrids.interpolate!(out, A, I.locator, I.geometry)
            @test maximum(abs, out .- A.data) < 1.0e-10
        end
    end
end

@testset "CubicInterpolator: partition of unity" begin
    # weights must sum to 1 everywhere, including the polar slots where the stencil degrades
    grid = OctahedralGaussianGrid(24)
    npoints = RingGrids.get_npoints(grid)
    londs, latds = RingGrids.get_londlatds(grid)

    I = CubicInterpolator(grid, npoints; NF = Float64)
    lon2 = mod.(londs .+ 3.7, 360)
    lat2 = clamp.(latds .+ 1.3, -90, 90)
    RingGrids.update_locator!(I, lon2, lat2)

    out = zeros(Float64, npoints)
    RingGrids.interpolate!(out, Field(ones(Float64, npoints), grid), I.locator, I.geometry)
    @test maximum(abs, out .- 1) < 1.0e-12

    # and explicitly at and beyond the outermost rings, where fewer than 4 rings are available
    I2 = CubicInterpolator(grid, 5; NF = Float64)
    RingGrids.update_locator!(I2, [0.0, 10.0, 20.0, 30.0, 40.0], [90.0, 89.9, -89.9, -90.0, 0.0])
    out2 = zeros(Float64, 5)
    RingGrids.interpolate!(out2, Field(ones(Float64, npoints), grid), I2.locator, I2.geometry)
    @test maximum(abs, out2 .- 1) < 1.0e-12

    # a smooth field must stay bounded at the poles, not extrapolate
    A = Field([sind(latds[i]) for i in 1:npoints], grid)
    RingGrids.interpolate!(out2, A, I2.locator, I2.geometry)
    @test all(-1.001 .<= out2 .<= 1.001)
end

@testset "CubicInterpolator: more accurate than AnvilInterpolator" begin
    grid = FullGaussianGrid(24)
    npoints = RingGrids.get_npoints(grid)
    londs, latds = RingGrids.get_londlatds(grid)

    # off-grid target points
    lon2 = mod.(londs .+ 3.7, 360)
    lat2 = clamp.(latds .+ 1.3, -90, 90)

    A = Field([cosd(latds[i])^4 * cosd(4londs[i]) for i in 1:npoints], grid)
    exact = [cosd(lat2[i])^4 * cosd(4lon2[i]) for i in 1:npoints]

    errors = Dict{Symbol, Float64}()
    for (name, IType) in ((:anvil, AnvilInterpolator), (:cubic, CubicInterpolator))
        I = IType(grid, npoints; NF = Float64)
        RingGrids.update_locator!(I, lon2, lat2)
        out = zeros(Float64, npoints)
        RingGrids.interpolate!(out, A, I.locator, I.geometry)
        errors[name] = sqrt(sum(abs2, out .- exact) / npoints)
    end

    # cubic should be at least an order of magnitude better on a smooth field
    @test errors[:cubic] < errors[:anvil] / 10
end

@testset "CubicInterpolator: damping under repeated shifts" begin
    # the semi-Lagrangian use case: the same field is interpolated every time step, so the
    # per-application amplitude loss compounds. This is what makes bilinear-class interpolation
    # unattractive for transport, and the main reason CubicInterpolator exists.
    grid = OctahedralGaussianGrid(48)
    npoints = RingGrids.get_npoints(grid)
    londs, latds = RingGrids.get_londlatds(grid)
    dlon = 360 / RingGrids.get_nlons(grid)[end] / 3      # a third of a grid cell

    amplitudes = Dict{Symbol, Float64}()
    for (name, IType) in ((:anvil, AnvilInterpolator), (:cubic, CubicInterpolator))
        I = IType(grid, npoints; NF = Float64)
        RingGrids.update_locator!(I, mod.(londs .+ dlon, 360), latds)

        A = Field([cosd(latds[i])^4 * cosd(4londs[i]) for i in 1:npoints], grid)
        amplitude0 = maximum(abs, A.data)
        buffer = zeros(Float64, npoints)
        for _ in 1:100
            RingGrids.interpolate!(buffer, A, I.locator, I.geometry)
            A.data .= buffer
        end
        amplitudes[name] = maximum(abs, A.data) / amplitude0
    end

    @test amplitudes[:cubic] > 0.99        # cubic retains essentially everything
    @test amplitudes[:anvil] < 0.95        # bilinear-class visibly damps
    @test amplitudes[:cubic] > amplitudes[:anvil]
end

@testset "Interpolation is interpolator-flexible" begin
    # the generic machinery must work for any AbstractInterpolator without special-casing:
    # `interpolator(grid, npoints; Interpolator)` picks the locator via `Locator(I)`, and
    # `update_locator!`/`interpolate!` dispatch on it.
    grid = OctahedralGaussianGrid(16)
    npoints = RingGrids.get_npoints(grid)
    londs, latds = RingGrids.get_londlatds(grid)

    for IType in (AnvilInterpolator, CubicInterpolator)
        I = RingGrids.interpolator(grid, npoints; Interpolator = IType)
        @test I isa IType
        @test I.locator isa RingGrids.Locator(IType)
        @test eltype(I) == RingGrids.DEFAULT_NF

        RingGrids.update_locator!(I, londs, latds)
        A = Field(rand(Float32, npoints), grid)
        out = zeros(Float32, npoints)
        RingGrids.interpolate!(out, A, I.locator, I.geometry)
        @test maximum(abs, out .- A.data) < 1.0f-5       # identity locator
    end

    # grid-to-grid interpolation with a non-default interpolator
    grid_out = FullGaussianGrid(12)
    I = RingGrids.interpolator(grid_out, grid; Interpolator = CubicInterpolator)
    A = Field([sind(latds[i]) for i in 1:npoints], grid)
    out = zeros(Float32, grid_out)
    RingGrids.interpolate!(out, A, I)
    _, latds_out = RingGrids.get_londlatds(grid_out)
    @test maximum(abs, out.data .- sind.(latds_out)) < 1.0f-3
end
