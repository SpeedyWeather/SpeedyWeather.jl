using NCDatasets

@testset "load_from_netcdf!" begin

    spectral_grid = SpectralGrid(trunc=31)
    (; GridVariable2D, GridVariable3D, grid) = spectral_grid

    @testset "input_data_path" begin
        p = SpeedyWeather.input_data_path("SpeedyWeather.jl/input_data", "orography.nc")
        @test isfile(p)

        p2 = SpeedyWeather.input_data_path("/tmp", "test.nc")
        @test p2 == "/tmp/test.nc"
    end

    @testset "2D: orography" begin
        # manual reference
        fullpath = SpeedyWeather.input_data_path("SpeedyWeather.jl/input_data", "orography.nc")
        ncfile = NCDataset(fullpath)
        orog_highres = FullGaussianField(ncfile["orog"].var[:, :], input_as=Matrix)
        orog_ref = zeros(GridVariable2D, grid)
        interp = RingGrids.interpolator(orog_ref, orog_highres, NF=Float32)
        RingGrids.interpolate!(orog_ref, orog_highres, interp)

        # via load_from_netcdf!
        orog_new = zeros(GridVariable2D, grid)
        load_from_netcdf!(orog_new, "SpeedyWeather.jl/input_data", "orography.nc", "orog")

        @test orog_ref ≈ orog_new
        @test maximum(abs, orog_new) > 0
    end

    @testset "3D: sea surface temperature" begin
        fullpath = SpeedyWeather.input_data_path("SpeedyWeather.jl/input_data", "sea_surface_temperature.nc")
        ncfile = NCDataset(fullpath)
        fill_value = ncfile["sst"].attrib["_FillValue"]
        sst_raw = FullGaussianGrid(ncfile["sst"].var[:, :, :], input_as=Matrix)
        sst_raw[sst_raw .=== fill_value] .= NaN
        sst_ref = zeros(GridVariable3D, grid, 12)
        interp = RingGrids.interpolator(sst_ref, sst_raw, NF=Float32)
        RingGrids.interpolate!(sst_ref, sst_raw, interp)

        sst_new = zeros(GridVariable3D, grid, 12)
        load_from_netcdf!(sst_new, "SpeedyWeather.jl/input_data", "sea_surface_temperature.nc", "sst";
            missing_value=NaN)

        mask = .!isnan.(sst_ref)
        @test sst_ref[mask] ≈ sst_new[mask]
        @test all(isnan.(sst_ref) .== isnan.(sst_new))
    end

    @testset "scale keyword" begin
        orog1 = zeros(GridVariable2D, grid)
        orog2 = zeros(GridVariable2D, grid)
        load_from_netcdf!(orog1, "SpeedyWeather.jl/input_data", "orography.nc", "orog")
        load_from_netcdf!(orog2, "SpeedyWeather.jl/input_data", "orography.nc", "orog"; scale=2)
        @test orog2 ≈ 2 .* orog1
    end
end
