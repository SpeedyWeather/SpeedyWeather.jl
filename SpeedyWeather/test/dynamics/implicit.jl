@testset "Implicit solver: S⁻¹ via LU decomposition and via Gauss-Jordan elimination (Reactant)" begin
    @testset for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, truncation = 32, nlayers = 8)
        model = PrimitiveDryModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, steps = 3)     # to have a non-trivial average temperature profile

        (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
        Tₖ = simulation.variables.dynamics.average_temperature_profile
        Δt = NF(2 * 1800 / model.planet.radius)    # some (scaled) time step

        implicit = deepcopy(model.implicit)
        SpeedyWeather.initialize!(implicit, Δt, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)
        implicit_gauss_jordan = deepcopy(implicit)
        ξ = implicit.centering * Δt
        SpeedyWeather.invert_S_gauss_jordan!(implicit_gauss_jordan, ξ)

        # S⁻¹ inverts S = 1 - ξ²λₗ(RL + UWᵀ) for every degree l, λₗ = -l(l-1) for 1-based l
        (; truncation, R, L, U, W) = implicit
        for l in 1:truncation
            S = one(R) - ξ^2 * (-l * (l - 1)) * (R * L + U * W')
            @test implicit.S⁻¹[l, :, :] * S ≈ one(S) rtol = 10 * sqrt(eps(NF))
        end

        # Gauss-Jordan elimination without pivoting vs LU decomposition
        @test implicit_gauss_jordan.S⁻¹[1:truncation, :, :] ≈ implicit.S⁻¹[1:truncation, :, :] rtol = 10 * sqrt(eps(NF))
    end
end
