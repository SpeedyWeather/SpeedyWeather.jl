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

@testset "Implicit operator L equals the linearised explicit temperature tendency" begin
    # For a resting, horizontally uniform atmosphere T = Tₖ at uniform surface pressure, the
    # explicit temperature tendency from a divergence field D is purely linear and has to equal
    # L*D, the linear terms that implicit_correction! moves between time steps
    @testset for T_top in (210, 290)    # non-isothermal and isothermal reference profile
        NF = Float64
        spectral_grid = SpectralGrid(; NF, truncation = 21, nlayers = 8)
        model = PrimitiveDryModel(spectral_grid; orography = NoOrography(spectral_grid), dynamics_only = true)
        model.feedback.verbose = false
        simulation = initialize!(model)
        vars = simulation.variables
        nlayers = spectral_grid.nlayers

        Tₖ = collect(range(NF(T_top), NF(290), nlayers))
        model.implicit.temp_profile .= Tₖ
        (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
        SpeedyWeather.initialize!(model.implicit, one(NF), Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)

        # zonal (m=0, real) divergence modes l=1:9, uniform T = Tₖ and pₛ = 1000 hPa, in both steps
        D = zeros(Complex{NF}, spectral_grid.spectrum, nlayers)
        D[2:10, :] .= 1.0e-6 .* randn(NF, 9, nlayers)
        for s in 1:2
            SpeedyWeather.get_step(vars.prognostic.vorticity, s) .= 0
            SpeedyWeather.get_step(vars.prognostic.divergence, s) .= D
            temp = SpeedyWeather.get_step(vars.prognostic.temperature, s)
            temp .= 0
            temp[1, :] .= Tₖ * sqrt(4π)     # constant c is c*sqrt(4π) in the l=m=0 harmonic
            pres = SpeedyWeather.get_step(vars.prognostic.pressure, s)
            pres .= 0
            pres[1] = log(NF(1.0e5)) * sqrt(4π)
        end

        SpeedyWeather.scale_prognostic!(vars, model.planet.radius)
        transform!(vars, model)
        transform!(vars, model)
        SpeedyWeather.reset_tendencies!(vars, model.time_stepping)
        SpeedyWeather.dynamics_tendencies!(vars, model)

        temp_tend = SpeedyWeather.get_tendency_step(vars.tendencies.temperature, model.time_stepping, SpeedyWeather.DynamicalCore())
        div = SpeedyWeather.get_step(vars.prognostic.divergence, 2)     # radius-scaled
        LD = transpose(model.implicit.L * transpose(div[2:10, :]))
        @test temp_tend[2:10, :] ≈ LD rtol = 1.0e-10
    end
end
