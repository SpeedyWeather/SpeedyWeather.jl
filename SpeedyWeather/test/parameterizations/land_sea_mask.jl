@testset "Land sea masks" begin
    @testset for Mask in (LandSeaMask, AquaPlanetMask, RockyPlanetMask)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
        mask = Mask(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; land_sea_mask = mask)
        simulation = initialize!(model)
        model.feedback.verbose = false
        run!(simulation, period = Day(5))
        @test simulation.model.feedback.nans_detected == false
    end
end

@testset "Land sea mask set! clamping" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    mask = LandSeaMask(spectral_grid)

    # out of [0, 1] range: should warn and clamp
    @test_logs (:warn,) set!(mask, 1.5)
    @test all(mask.land_fraction .== 1)

    # already in [0, 1] range: no warning, values unchanged
    set!(mask, 0.3)
    @test_logs set!(mask, 0.3)
    @test all(mask.land_fraction .== eltype(mask.land_fraction)(0.3))
end
