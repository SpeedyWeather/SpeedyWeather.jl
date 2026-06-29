# `_vertical_advection!(::GPU, ...)` dispatches between two kernels depending on npoints
# (vertical_advection_kernel!, one thread per (ij,k); vertical_advection_column_kernel!, one
# thread per ij looping over k), checked here directly against real GPU data for every
# scheme, bypassing the npoints threshold so both kernels are exercised regardless of this
# test grid's resolution. Uses ≈, not ===: under --check-bounds=yes the inserted bounds
# checks change the generated code enough to shift the GPU compiler's FP reordering/FMA
# decisions between the two differently-shaped kernels by a few ULPs for Upwind/WENO
# (max relative diff ~2e-4, well within tolerance) -- the same class of GPU FP
# non-determinism vertical_integration.jl's GPU test already accounts for with ≈.
@testset "Vertical advection GPU kernels" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(; architecture = arch)

    advection_schemes = (
        SpeedyWeather.CenteredVerticalAdvection(spectral_grid),
        SpeedyWeather.UpwindVerticalAdvection(spectral_grid),
        SpeedyWeather.WENOVerticalAdvection(spectral_grid),
    )

    for advection_scheme in advection_schemes
        model = PrimitiveWetModel(spectral_grid; vertical_advection = advection_scheme, dynamics_only = true)
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 2)

        vars, model = SpeedyWeather.unpack(simulation)
        Δσ = model.geometry.σ_levels_thick
        w = vars.dynamics.w
        ξ = vars.grid.u
        nlayers = size(ξ, 2)
        npoints = size(ξ, 1)

        ξ_tend_pointwise = deepcopy(vars.tendencies.grid.u)
        ξ_tend_column = deepcopy(vars.tendencies.grid.u)

        SpeedyWeather.launch!(
            arch, SpeedyWeather.RingGridWorkOrder, (npoints, nlayers),
            SpeedyWeather.vertical_advection_kernel!,
            ξ_tend_pointwise, 1, w, ξ, 1, Δσ, nlayers, advection_scheme
        )
        SpeedyWeather.launch!(
            arch, SpeedyWeather.LinearWorkOrder, (npoints,),
            SpeedyWeather.vertical_advection_column_kernel!,
            ξ_tend_column, 1, w, ξ, 1, Δσ, nlayers, advection_scheme
        )

        @test ξ_tend_pointwise.data ≈ ξ_tend_column.data
    end
end
