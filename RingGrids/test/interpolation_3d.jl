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
