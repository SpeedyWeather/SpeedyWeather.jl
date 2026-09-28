# minimal duck-typed stand-in for SpeedyWeather.Particle (not a RingGrids dependency)
# interpolate_3D! only ever accesses the .σ field of `positions`
struct TestParticle3D{NF}
    σ::NF
end

@testset "3D interpolation: vertical profile" begin
    npoints = 50
    nlayers = 6

    @testset for Grid in (
            FullGaussianGrid,
            OctahedralGaussianGrid,
            HEALPixGrid,
        )

        @testset for NF in (Float32, Float64)
            grid = Grid(8)
            σ_levels_full = collect(range(NF(0.1), NF(0.9), length = nlayers))

            # field constant across the horizontal, linear in σ across layers
            A = zeros(NF, grid, nlayers)
            for k in 1:nlayers
                A[:, k] .= σ_levels_full[k]
            end

            geometry = RingGrids.GridGeometry(A)
            locator = RingGrids.AnvilLocator(NF, npoints, nlayers)

            λs = 360 * rand(NF, npoints)
            θs = 180 * rand(NF, npoints) .- 90
            RingGrids.update_locator!(locator, geometry, λs, θs)

            # random σ within [σ_levels_full[1], σ_levels_full[end]] so no pinning occurs
            σs = σ_levels_full[1] .+ (σ_levels_full[end] - σ_levels_full[1]) .* rand(NF, npoints)
            positions = [TestParticle3D(σ) for σ in σs]

            Aout = zeros(NF, npoints)
            RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels_full)

            # A is horizontally constant and piecewise-linear in σ at σ_levels_full,
            # so the vertically-blended interpolation should reconstruct σ itself
            @test Aout ≈ σs
        end
    end
end

@testset "3D interpolation: pin outside σ range" begin
    NF = Float32
    grid = FullGaussianGrid(8)
    nlayers = 4
    σ_levels_full = NF[0.1, 0.3, 0.6, 0.9]

    A = zeros(NF, grid, nlayers)
    for k in 1:nlayers
        A[:, k] .= σ_levels_full[k]
    end

    geometry = RingGrids.GridGeometry(A)
    locator = RingGrids.AnvilLocator(NF, 2, nlayers)
    RingGrids.update_locator!(locator, geometry, NF[10, 200], NF[0, 0])

    # below and above the σ range: should be pinned, not extrapolated
    positions = [TestParticle3D(NF(-1)), TestParticle3D(NF(2))]
    Aout = zeros(NF, 2)
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels_full)

    @test Aout[1] ≈ σ_levels_full[1]
    @test Aout[2] ≈ σ_levels_full[end]
end

@testset "3D interpolation: pole averaging" begin
    NF = Float32
    grid = FullGaussianGrid(8)
    nlayers = 3
    σ_levels_full = NF[0.2, 0.5, 0.8]

    # field = ring latitude in degrees, identical on every layer
    A = zeros(NF, grid, nlayers)
    geometry = RingGrids.GridGeometry(A)
    for (j, ring) in enumerate(RingGrids.eachring(grid))
        θ = geometry.latd[j + 1]
        for ij in ring
            A[ij, :] .= θ
        end
    end

    locator = RingGrids.AnvilLocator(NF, 1, nlayers)
    RingGrids.update_locator!(locator, geometry, NF[0], NF[89.99])   # just south of the north pole

    # exactly on the 2nd σ level so no vertical blending occurs either
    positions = [TestParticle3D(σ_levels_full[2])]
    Aout = zeros(NF, 1)
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels_full)

    # north pole value is the average of the first ring, i.e. its (constant) latitude
    @test Aout[1] ≈ geometry.latd[2]
end

@testset "3D interpolation: size checks" begin
    NF = Float32
    grid = FullGaussianGrid(8)
    nlayers = 4
    A = zeros(NF, grid, nlayers)
    geometry = RingGrids.GridGeometry(A)
    σ_levels_full = NF[0.1, 0.3, 0.6, 0.9]
    positions = [TestParticle3D(NF(0.5)) for _ in 1:3]
    Aout = zeros(NF, 3)

    # default locator has pole buffers for a single layer only
    locator_single_layer = RingGrids.AnvilLocator(NF, 3)
    @test_throws DimensionMismatch RingGrids.interpolate_3D!(Aout, A, locator_single_layer, geometry, positions, σ_levels_full)

    locator = RingGrids.AnvilLocator(NF, 3, nlayers)
    @test_throws DimensionMismatch RingGrids.interpolate_3D!(zeros(NF, 2), A, locator, geometry, positions, σ_levels_full)
    @test_throws DimensionMismatch RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions[1:2], σ_levels_full)

    # faces need nlayers + 1 σ levels, full levels aren't enough
    @test_throws DimensionMismatch RingGrids.interpolate_3D!(
        Aout, A, locator, geometry, positions, SigmaFaceBelow(σ_levels_full, zero(NF))
    )
end

@testset "3D interpolation: single layer" begin
    NF = Float32
    grid = FullGaussianGrid(8)
    A = zeros(NF, grid, 1)
    A[:, 1] .= 3
    geometry = RingGrids.GridGeometry(A)
    locator = RingGrids.AnvilLocator(NF, 2, 1)
    RingGrids.update_locator!(locator, geometry, NF[0, 180], NF[45, -45])

    # a single level means the field is constant in the vertical, for any σ
    positions = [TestParticle3D(NF(0.2)), TestParticle3D(NF(0.9))]
    Aout = zeros(NF, 2)
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, NF[0.5])
    @test Aout ≈ [3, 3]
end
