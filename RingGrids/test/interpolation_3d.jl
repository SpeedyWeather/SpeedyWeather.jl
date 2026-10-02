# minimal duck-typed stand-in for SpeedyWeather.Particle (not a RingGrids dependency)
# interpolate_3D! only ever accesses the .σ field of `positions`
struct TestParticle3D{NF}
    σ::NF
end

@testset "find_vertical_bracket" begin
    σf = Float32[0.1, 0.3, 0.5, 0.7, 0.9]

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(0.0f0, σf)
    @test k_lo == 1 && k_hi == 2 && α == 0.0f0       # below grid top → pin

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(σf[1], σf)
    @test k_lo == 1 && k_hi == 2 && α == 0.0f0       # exactly at top layer

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(0.4f0, σf)
    @test k_lo == 2 && k_hi == 3 && α ≈ 0.5f0        # midpoint between layers 2 and 3

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(σf[end], σf)
    @test k_lo == 4 && k_hi == 5 && α == 1.0f0       # exactly at bottom layer

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(1.0f0, σf)
    @test k_lo == 4 && k_hi == 5 && α == 1.0f0       # below grid bottom → pin

    @test all(
        0 ≤ α ≤ 1 for (_, _, α) in
            [RingGrids.find_vertical_bracket(σ, σf) for σ in range(0.0f0, 1.0f0, 50)]
    )

    # single level: that level with zero weight, for any σ
    @test RingGrids.find_vertical_bracket(0.2f0, Float32[0.5]) == (1, 1, 0.0f0)
    @test RingGrids.find_vertical_bracket(0.9f0, Float32[0.5]) == (1, 1, 0.0f0)
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
        Aout, A, locator, geometry, positions, SigmaFaceBelow(σ_levels_full)
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

@testset "3D interpolation: face staggering" begin
    NF = Float32
    grid = FullGaussianGrid(8)
    nlayers = 4
    σ_half = NF[0, 0.2, 0.5, 0.7, 1]
    f(σ) = 1 + σ * (1 - σ)                      # equal to the boundary value 1 at σ = 0 and 1

    # expected: piecewise linear interpolation of f between half levels
    function expected(σ)
        k = findlast(≤(σ), σ_half[1:(end - 1)])
        α = (σ - σ_half[k]) / (σ_half[k + 1] - σ_half[k])
        return f(σ_half[k]) + (f(σ_half[k + 1]) - f(σ_half[k])) * α
    end

    geometry = RingGrids.GridGeometry(zeros(NF, grid, nlayers))
    locator = RingGrids.AnvilLocator(NF, 6, nlayers)
    RingGrids.update_locator!(locator, geometry, NF[0, 60, 120, 180, 240, 300], NF[80, 45, 10, -10, -45, -80])
    σs = NF[0, 0.1, 0.2, 0.35, 0.95, 1]
    positions = [TestParticle3D(σ) for σ in σs]

    # FaceBelow: layer k stores σ_half[k+1], top (σ=0) not stored
    A = zeros(NF, grid, nlayers)
    for k in 1:nlayers
        A[:, k] .= f(σ_half[k + 1])
    end
    Aout = zeros(NF, 6)
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, SigmaFaceBelow(σ_half, one(NF)))
    @test Aout ≈ expected.(σs)

    # FaceAbove: layer k stores σ_half[k], bottom (σ=1) not stored
    for k in 1:nlayers
        A[:, k] .= f(σ_half[k])
    end
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, positions, SigmaFaceAbove(σ_half, one(NF)))
    @test Aout ≈ expected.(σs)

    # one-argument constructors default to a zero boundary value of the σ levels' number type
    @test SigmaFaceBelow(σ_half).top_boundary_condition === zero(NF)
    @test SigmaFaceAbove(σ_half).bottom_boundary_condition === zero(NF)
end

@testset "interpolate! forwards to interpolate_3D!" begin
    NF = Float32
    nlayers = 4
    npoints = 20
    A = rand(NF, FullGaussianGrid(8), nlayers)
    geometry = RingGrids.GridGeometry(A)
    locator = RingGrids.AnvilLocator(NF, npoints, nlayers)
    RingGrids.update_locator!(locator, geometry, 360 * rand(NF, npoints), 180 * rand(NF, npoints) .- 90)
    positions = [TestParticle3D(σ) for σ in rand(NF, npoints)]

    σ_full = NF[0.1, 0.3, 0.6, 0.9]
    σ_half = NF[0, 0.2, 0.5, 0.7, 1]
    @testset for vertical in (SigmaCenter(σ_full), σ_full, SigmaFaceBelow(σ_half))
        out1, out2 = zeros(NF, npoints), zeros(NF, npoints)
        interpolate!(out1, A, locator, geometry, positions, vertical)
        interpolate_3D!(out2, A, locator, geometry, positions, vertical)
        @test out1 == out2
    end
end
