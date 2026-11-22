# The safety checks on the latitudes in `update_locator!` can't be traced, so with Reactant the
# locator is always updated without them (`unsafe = true`, RingGridsReactantExt) and with a
# branchless kernel instead of the binary search. Same locator as on the CPU, eagerly and compiled
# (as in the particle advection), for valid latitudes and invalid ones (clamped to the poles)

@testset "update_locator! on Reactant" begin
    RingGrids = SpeedyWeather.RingGrids
    arch = ReactantDevice()
    NF = Float32
    locator_ints(L) = (L.js, L.ij_as, L.ij_bs, L.ij_cs, L.ij_ds)
    locator_floats(L) = (L.Δys, L.Δabs, L.Δcds)

    function test_same_locator(L, L_cpu)
        @test all(map((a, b) -> Array(a) == b, locator_ints(L), locator_ints(L_cpu)))
        @test all(map((a, b) -> Array(a) ≈ b, locator_floats(L), locator_floats(L_cpu)))
        return nothing
    end

    for Grid in (OctahedralGaussianGrid, FullGaussianGrid, HEALPixGrid)
        grid_cpu = Grid(8)
        grid = Grid(8, arch)

        # valid latitudes incl the poles and on the rings, then invalid ones
        latd = NF.(RingGrids.get_latd(grid_cpu))
        θs = vcat(NF[90, 0, -90], latd, 180 * rand(NF, 50) .- 90, NF[90.5, -90.5, 1000, -1000])
        λs = 720 * rand(NF, length(θs)) .- 360

        I_cpu = RingGrids.interpolator(grid_cpu, length(θs); NF)
        RingGrids.update_locator!(I_cpu, λs, θs, unsafe = true)
        @test I_cpu.locator.js[(end - 3):end] == [0, I_cpu.geometry.nlat, 0, I_cpu.geometry.nlat]

        λs_device, θs_device = on_architecture(arch, λs), on_architecture(arch, θs)

        # eagerly, e.g. creating an interpolator outside of compiled code
        I = RingGrids.interpolator(grid, length(θs); NF)
        RingGrids.update_locator!(I, λs_device, θs_device)
        test_same_locator(I.locator, I_cpu.locator)

        # compiled, e.g. the particle advection inside the time step
        I = RingGrids.interpolator(grid, length(θs); NF)
        Reactant.@jit RingGrids.update_locator!(I.locator, I.geometry, λs_device, θs_device)
        test_same_locator(I.locator, I_cpu.locator)
    end
end
