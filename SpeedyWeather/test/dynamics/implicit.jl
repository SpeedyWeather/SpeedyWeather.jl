@testset "Traceable initialization of the implicit solver matches initialize!" begin
    @testset for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, truncation = 32, nlayers = 8)
        model = PrimitiveDryModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, steps = 3)     # to have a non-trivial average temperature profile

        (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
        Tₖ = simulation.variables.dynamics.average_temperature_profile
        Δt = NF(2 * 1800 / model.planet.radius)    # some (scaled) time step

        implicit = deepcopy(model.implicit)
        implicit_traceable = deepcopy(model.implicit)
        SpeedyWeather.initialize!(implicit, Δt, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)
        SpeedyWeather.initialize_traceable!(implicit_traceable, Δt, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)

        for name in (:temp_profile, :R, :U, :L, :W)
            @test getfield(implicit_traceable, name) ≈ getfield(implicit, name)
        end

        # Gauss-Jordan elimination without pivoting vs LU decomposition
        (; truncation) = implicit
        @test implicit_traceable.S⁻¹[1:truncation, :, :] ≈ implicit.S⁻¹[1:truncation, :, :] rtol = 10 * sqrt(eps(NF))
    end
end
