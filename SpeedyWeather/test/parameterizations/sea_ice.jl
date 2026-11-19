@testset "Sea ice thickness, volume and concentration" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    sea_ice = ThermodynamicSeaIce(spectral_grid)
    (; thin_ice_thickness, pack_ice_thickness) = sea_ice

    @test SpeedyWeather.sea_ice_thickness(0, sea_ice) == thin_ice_thickness
    @test SpeedyWeather.sea_ice_thickness(1, sea_ice) == pack_ice_thickness

    # the concentration is the exact inverse of the volume, also for a fixed thickness (pack = thin)
    @testset for sea_ice in (sea_ice, ThermodynamicSeaIce(spectral_grid, pack_ice_thickness = 0.2))
        ℵs = Float32.(0:0.01:1)
        Vs = SpeedyWeather.sea_ice_volume.(ℵs, Ref(sea_ice))
        @test all(diff(Vs) .> 0)
        @test SpeedyWeather.sea_ice_concentration.(Vs, Ref(sea_ice)) ≈ ℵs atol = 1.0e-6
        @test SpeedyWeather.sea_ice_concentration(-1.0f0, sea_ice) == 0
        @test SpeedyWeather.sea_ice_concentration(10.0f0, sea_ice) == 1
    end

    # beyond full cover the ice thickens: h = V, continuous at V = pack_ice_thickness
    h(V) = SpeedyWeather.sea_ice_thickness_from_volume(V, sea_ice)
    @test h(4.0f0) ≈ 4
    @test h(pack_ice_thickness) ≈ pack_ice_thickness
    @test h(pack_ice_thickness - 1.0f-4) ≈ pack_ice_thickness atol = 1.0e-3
    @test h(0.0f0) == thin_ice_thickness
end

@testset "Sea ice surface energy balance" begin
    ϵσ = 0.98f0 * 5.67f-8
    T_f, T_m = 271.35f0, 273.15f0

    # Newton against bisection, from thin (stiff, G = 10 W/m²/K) to thick ice, cold and mild forcing
    @testset for G in (10.0f0, 2 / 2.5f0, 0.2f0), A in (150.0f0, 250.0f0)
        T_i, F_c, surplus = SpeedyWeather.sea_ice_surface_balance(A, G, ϵσ, T_f, T_m, 230.0f0)
        f(T) = ϵσ * T^4 + G * (T - T_f) - A
        lo, hi = 100.0, 400.0
        for _ in 1:60
            mid = (lo + hi) / 2
            f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        @test T_i ≈ min(lo, T_m) atol = 1.0e-3
        @test F_c ≈ G * (T_f - T_i)
        @test surplus ≈ max(A + F_c - ϵσ * T_i^4, 0) atol = 1.0e-3
        @test surplus >= 0                          # never negative, zero unless capped
    end

    # warm forcing: capped at melting, the surplus melts ice
    T_i, F_c, surplus = SpeedyWeather.sea_ice_surface_balance(500.0f0, 0.8f0, ϵσ, T_f, T_m, 250.0f0)
    @test T_i == T_m
    @test F_c < 0                                   # conducted down into the ocean from the melting surface
    @test surplus > 0

    # very conductive ice: the surface is at the freezing temperature, no conduction with T_i = T_f
    T_i, F_c, _ = SpeedyWeather.sea_ice_surface_balance(200.0f0, 1.0f6, ϵσ, T_f, T_m, 250.0f0)
    @test T_i ≈ T_f atol = 1.0e-2

    # Float32 in, Float32 out
    @test all(x -> x isa Float32, SpeedyWeather.sea_ice_surface_balance(200.0f0, 1.0f0, ϵσ, T_f, T_m, 250.0f0))
end

@testset "Sea ice freezing and melting conserve energy" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    model = PrimitiveWetModel(spectral_grid; ocean = SlabOcean(spectral_grid))
    simulation = initialize!(model)
    vars = simulation.variables
    sea_ice = model.sea_ice
    sst = vars.prognostic.ocean.sea_surface_temperature
    V = vars.prognostic.ocean.sea_ice_volume
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    ocean = model.land_sea_mask.land_fraction.data .< 1
    C = model.ocean.heat_capacity_mixed_layer
    ρL = sea_ice.density * sea_ice.latent_heat_fusion
    T_f = sea_ice.freezing_temperature
    energy() = Float64.(C .* sst.data .- ρL .* V.data)[ocean, :]   # ocean heat minus latent heat in ice

    # freezing from open water, 1 K below freezing: as main's freeze rate (0.5 per K) at first freeze-up
    sst.data .= T_f - 1
    V.data .= 0
    E0 = energy()
    filter!(vars, sea_ice, model)
    @test energy() ≈ E0 rtol = 1.0e-6
    @test all(sst.data[ocean, :] .== T_f)
    @test all(x -> isapprox(x, 0.5, atol = 0.01), ℵ.data[ocean, :])

    # freezing beyond full cover thickens the ice, the water never stays below freezing
    sst.data .= T_f - 5
    V.data .= sea_ice.pack_ice_thickness
    filter!(vars, sea_ice, model)
    @test all(V.data[ocean, :] .> sea_ice.pack_ice_thickness)
    @test all(ℵ.data[ocean, :] .== 1)
    @test all(sst.data[ocean, :] .== T_f)

    # bottom melt: the ocean pays the latent heat, never more ice than present
    sst.data .= T_f + 2
    V.data .= 0.1
    E0 = energy()
    filter!(vars, sea_ice, model)
    @test energy() ≈ E0 rtol = 1.0e-6
    @test all(0 .<= V.data[ocean, :] .< 0.1)
    @test all(sst.data[ocean, :] .< T_f + 2)

    # the basal flux γ (SST - T_f) melts exactly its energy in ice volume
    Δt = SpeedyWeather.time_step(model.time_stepping, :ocean, vars.prognostic.clock)
    sst.data .= T_f + 1
    V.data .= 1
    filter!(vars, sea_ice, model)
    @test 1 .- V.data[ocean, :] ≈ fill(sea_ice.basal_melt_coefficient * Δt / ρL, count(ocean), size(V, 2)) rtol = 1.0e-3

    # as ℵ → 0 (dℵ/dV = 1/thin_ice_thickness) that is main's melt rate of 1e-7 per K per second
    @test sea_ice.basal_melt_coefficient / (ρL * sea_ice.thin_ice_thickness) ≈ 1.0e-7 rtol = 0.01

    # top melt beyond the ice present (negative volume after the time step) warms the ocean instead
    sst.data .= T_f + 1
    V.data .= -0.01
    E0 = energy()
    filter!(vars, sea_ice, model)
    @test energy() ≈ E0 rtol = 1.0e-6
    @test all(V.data[ocean, :] .== 0)
    @test all(ℵ.data[ocean, :] .== 0)

    # no change above freezing without ice, nor at freezing with ice
    sst.data .= T_f + 1
    V.data .= 0
    filter!(vars, sea_ice, model)
    @test all(sst.data[ocean, :] .== T_f + 1)
    sst.data .= T_f
    V.data .= 1
    filter!(vars, sea_ice, model)
    @test all(V.data[ocean, :] .== 1)

    # over a prescribed ocean the ice grows with a nominal 50 m mixed layer, the SST is unchanged
    model = PrimitiveWetModel(spectral_grid; ocean = SeasonalOceanClimatology(spectral_grid))
    vars = initialize!(model).variables
    @test !haskey(vars.tendencies.ocean, :sea_surface_temperature)
    sst = vars.prognostic.ocean.sea_surface_temperature
    vars.prognostic.ocean.sea_ice_volume .= 0
    sst.data .= T_f - 1
    filter!(vars, model.sea_ice, model)
    ocean = model.land_sea_mask.land_fraction.data .< 1
    @test all(sst.data[ocean, :] .== T_f - 1)
    @test vars.prognostic.ocean.sea_ice_volume.data[ocean, :] ≈
        fill(model.sea_ice.prescribed_ocean_heat_capacity / ρL, count(ocean), 1) rtol = 1.0e-5

    # the prescribed SST stays below freezing: the ice grows to full cover but never thicker (bounded, as main's clamp)
    for _ in 1:100
        filter!(vars, model.sea_ice, model)
    end
    @test all(vars.prognostic.ocean.sea_ice_volume.data[ocean, :] .≈ model.sea_ice.pack_ice_thickness)
    @test all(vars.prognostic.ocean.sea_ice_concentration.data[ocean, :] .≈ 1)
    # thicker ice (e.g. from an earlier interactive run) is kept, not thinned
    vars.prognostic.ocean.sea_ice_volume.data[ocean, :] .= 4
    filter!(vars, model.sea_ice, model)
    @test all(vars.prognostic.ocean.sea_ice_volume.data[ocean, :] .== 4)
end

@testset "Thickness helpers with pack thinner than thin ice" begin
    # no NaN if pack_ice_thickness < thin_ice_thickness: pack is taken as at least thin
    for V in (0.0f0, 0.1f0, 0.5f0, 3.0f0)
        ℵ = SpeedyWeather.sea_ice_concentration(V, 0.5f0, 0.2f0)
        @test isfinite(ℵ) && 0 <= ℵ <= 1
        @test isfinite(SpeedyWeather.sea_ice_thickness_from_volume(V, 0.5f0, 0.2f0))
    end
end

@testset "Slab ocean pays the top melt" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    model = PrimitiveWetModel(spectral_grid; ocean = SlabOcean(spectral_grid))
    vars = initialize!(model).variables
    ocean = model.land_sea_mask.land_fraction.data .< 1
    dsst = vars.tendencies.ocean.sea_surface_temperature

    vars.parameterizations.ocean.sea_ice_top_melt .= 0
    SpeedyWeather.timestep!(vars, model.ocean, model)
    dsst0 = copy(dsst.data)

    vars.parameterizations.ocean.sea_ice_top_melt .= 10
    SpeedyWeather.timestep!(vars, model.ocean, model)
    @test (dsst0 .- dsst.data)[ocean] ≈ fill(10 / model.ocean.heat_capacity_mixed_layer, count(ocean)) rtol = 1.0e-3

    # and the sea ice volume tendency is the top melt in ice volume
    SpeedyWeather.timestep!(vars, model.sea_ice, model)
    dV = vars.tendencies.ocean.sea_ice_volume.data
    @test dV[ocean, 1] ≈ fill(-10 / (model.sea_ice.density * model.sea_ice.latent_heat_fusion), count(ocean)) rtol = 1.0e-5
end

@testset "Downward longwave needed in advance" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    # only ThermodynamicSeaIce solves an ice surface balance that needs the surface downward longwave first
    @test SpeedyWeather.needs_surface_longwave_down(ThermodynamicSeaIce(spectral_grid))
    @test !SpeedyWeather.needs_surface_longwave_down(PrescribedSeaIce())
    @test !SpeedyWeather.needs_surface_longwave_down(nothing)

    # the pre-pass gives bitwise the surface downward longwave of OneBandLongwave's own downward beam
    model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave = OneBandLongwave(spectral_grid)))
    vars = initialize!(model).variables
    longwave = model.radiation.longwave
    vars.parameterizations.surface_pressure .= 1.0e5
    for k in 1:spectral_grid.nlayers
        vars.grid.temperature[:, k, :] .= 210 + 15 * (k - 1)
    end
    vars.grid.humidity .= 0.004
    ij = 1
    SpeedyWeather.parameterization!(ij, vars, longwave, model)
    T = SpeedyWeather.get_prognostic_step(vars.grid.temperature, model.time_stepping, longwave.radiative_transfer)
    t = SpeedyWeather.transmissivity!(ij, vars, longwave.transmissivity, model)     # same column, same values
    σ = model.atmosphere.stefan_boltzmann
    D = SpeedyWeather.surface_longwave_down(ij, T, t, σ, spectral_grid.nlayers)
    @test D == vars.parameterizations.surface_longwave_down[ij]
end

@testset "Sea ice output variables" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    vars = simulation.variables
    vars.prognostic.ocean.sea_ice_volume .= 1.5
    vars.parameterizations.ocean.sea_ice_surface_temperature .= 250
    volume = OutputVariables.SeaIceVolume()
    surface_temperature = OutputVariables.SeaIceSurfaceTemperature()
    @test volume isa SpeedyWeather.AbstractOutputVariable
    @test OutputVariables.path(volume, simulation) === vars.prognostic.ocean.sea_ice_volume
    @test OutputVariables.path(surface_temperature, simulation) === vars.parameterizations.ocean.sea_ice_surface_temperature
    @test surface_temperature.transform(250.0f0) ≈ 250 - 273.15     # written in ˚C
    @test volume.unit == "m"
end
