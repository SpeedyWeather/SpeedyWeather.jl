@testset "Saturation vapor pressure formulations" begin
    atmosphere = EarthAtmosphere(Float64)

    # reference values [Pa] (Murphy and Koop, 2005), over liquid water
    liquid = [(253.15, 125.6), (273.15, 611.2), (293.15, 2339.0), (303.15, 4246.0), (313.15, 7385.0)]

    # over ice
    ice = [(233.15, 12.84), (243.15, 37.98)]

    # accurate formulations over liquid within 0.5% (IFS within 1% at -20˚C as mixed phase)
    for SVP in (ClausiusClapeyronLinearL, MagnusBolton)
        svp = SVP(Float64)
        for (T, e) in liquid
            @test SpeedyWeather.saturation_vapor_pressure(svp, T, atmosphere) ≈ e rtol = 5.0e-3
        end
    end

    svp = IFSMixedPhase(Float64)
    for (T, e) in liquid[2:end]
        @test SpeedyWeather.saturation_vapor_pressure(svp, T, atmosphere) ≈ e rtol = 5.0e-3
    end
    for (T, e) in ice
        @test SpeedyWeather.saturation_vapor_pressure(svp, T, atmosphere) ≈ e rtol = 1.0e-2
    end

    # constant L overestimates at warm temperatures (2-5% at 30-40˚C)
    svp = ClausiusClapeyronConstantL(Float64)
    @test SpeedyWeather.saturation_vapor_pressure(svp, 273.15, atmosphere) == 610.78
    @test 1.01 < SpeedyWeather.saturation_vapor_pressure(svp, 303.15, atmosphere) / 4246 < 1.05
end

@testset "Saturation vapor pressure gradients" begin
    atmosphere = EarthAtmosphere(Float64)

    # avoid the IFS phase boundaries 250.16K and 273.16K where the derivative is only continuous in α
    temperatures = 213.15:2.5:318.15
    h = 1.0e-4
    for SVP in (ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton)
        svp = SVP(Float64)
        for T in temperatures
            e⁺ = SpeedyWeather.saturation_vapor_pressure(svp, T + h, atmosphere)
            e⁻ = SpeedyWeather.saturation_vapor_pressure(svp, T - h, atmosphere)
            @test SpeedyWeather.grad(svp, T, atmosphere) ≈ (e⁺ - e⁻) / 2h rtol = 1.0e-6
        end
    end

    # IFS mixed phase is continuous at both ends of the mixed phase
    svp = IFSMixedPhase(Float64)
    for T in (svp.T_ice, svp.T₀)
        e⁺ = SpeedyWeather.saturation_vapor_pressure(svp, T + 1.0e-9, atmosphere)
        e⁻ = SpeedyWeather.saturation_vapor_pressure(svp, T - 1.0e-9, atmosphere)
        @test e⁺ ≈ e⁻ rtol = 1.0e-6
    end
end

@testset "Saturation vapor pressure in EarthAtmosphere" begin
    spectral_grid = SpectralGrid(truncation = 16, nlayers = 2)

    # default is the IFS mixed phase formulation
    @test EarthAtmosphere(spectral_grid).saturation isa IFSMixedPhase{Float32}

    # Clausius-Clapeyron with constant L is bit-identical to the previous default
    atmosphere = EarthAtmosphere(spectral_grid, saturation = ClausiusClapeyronConstantL(spectral_grid))
    for T in (230.0f0, 273.15f0, 300.0f0)
        e_old = 610.78f0 * exp(atmosphere.latent_heat_condensation / atmosphere.R_vapor * (inv(273.15f0) - inv(T)))
        @test SpeedyWeather.saturation_vapor_pressure(T, atmosphere) === e_old
        @test SpeedyWeather.saturation_humidity(T, 1.0f5, atmosphere) ===
            atmosphere.mol_ratio * e_old / 1.0f5
    end

    # swap formulation, Float32 in and out
    for SVP in (ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton)
        atmosphere = EarthAtmosphere(spectral_grid, saturation = SVP(spectral_grid))
        @test atmosphere.saturation isa SVP{Float32}
        @test SpeedyWeather.saturation_vapor_pressure(260.0f0, atmosphere) isa Float32
        @test SpeedyWeather.grad_saturation_vapor_pressure(260.0f0, atmosphere) isa Float32
        @test SpeedyWeather.grad_saturation_vapor_pressure(260.0f0, atmosphere) > 0
    end

    # saturation humidity gradient consistent with the formulation, in Float64 vs finite differences
    for SVP in (ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton)
        atmosphere = EarthAtmosphere(Float64, saturation = SVP(Float64))
        for T in (240.0, 265.0, 290.0, 305.0), p in (3.0e4, 1.0e5)
            h = 1.0e-4
            q⁺ = SpeedyWeather.saturation_humidity(T + h, p, atmosphere)
            q⁻ = SpeedyWeather.saturation_humidity(T - h, p, atmosphere)
            @test SpeedyWeather.grad_saturation_humidity(T, p, atmosphere) ≈ (q⁺ - q⁻) / 2h rtol = 1.0e-6
        end
    end

    # for constant L this is the previously hardcoded Clausius-Clapeyron derivative
    atmosphere = EarthAtmosphere(Float64, saturation = ClausiusClapeyronConstantL(Float64))
    (; latent_heat_condensation, R_vapor) = atmosphere
    for T in (250.0, 280.0, 300.0)
        qₛ = SpeedyWeather.saturation_humidity(T, 8.0e4, atmosphere)
        @test SpeedyWeather.grad_saturation_humidity(T, 8.0e4, atmosphere) ≈
            qₛ * latent_heat_condensation / (R_vapor * T^2)
    end

    # also through the model constructor
    atmosphere = EarthAtmosphere(spectral_grid, saturation = MagnusBolton(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; atmosphere)
    @test model.atmosphere.saturation isa MagnusBolton{Float32}

    # dry atmosphere has no saturation
    dry = EarthDryAtmosphere(spectral_grid)
    @test SpeedyWeather.saturation_vapor_pressure(280.0f0, dry) == 0
    @test SpeedyWeather.grad_saturation_vapor_pressure(280.0f0, dry) == 0
    @test SpeedyWeather.grad_saturation_humidity(280.0f0, 1.0f5, dry) == 0
end
