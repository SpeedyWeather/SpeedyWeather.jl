using Terrarium
using Test
using Dates
using NCDatasets
using Zarr

# The SpeedyWeatherTerrariumExt extension activates once both packages are loaded.
const SWTerrariumOutputExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherTerrariumExt)
@assert SWTerrariumOutputExt !== nothing "SpeedyWeatherTerrariumExt failed to load"
const TerrariumOutputVariable = SWTerrariumOutputExt.TerrariumOutputVariable

@testset "TerrariumOutput" begin
    # Small coupled setup as in test/parameterizations/terrarium_coupling.jl
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    Nz = 4
    land_mask = land_sea_mask.land_fraction .> 0
    column_grid = Terrarium.ColumnRingGrid(
        Terrarium.CPU(), Float32,
        Terrarium.ExponentialSpacing(; N = Nz, Δz_min = 0.05),
        ring_grid,
        land_mask,
    )

    soil = Terrarium.SoilEnergyWaterCarbon(
        eltype(column_grid);
        hydrology = Terrarium.SoilHydrology(eltype(column_grid)),
    )
    terrarium_model = Terrarium.LandModel(
        column_grid;
        initializer = Terrarium.SoilInitializer(eltype(column_grid)),
        vegetation = nothing,
        soil,
    )

    @testset "constructors + metadata" begin
        # 3D (subsurface) variable: metadata derived from Terrarium descriptors
        st = TerrariumOutput(terrarium_model, :temperature)
        @test st isa TerrariumOutputVariable
        @test st isa SpeedyWeather.AbstractOutputVariable
        @test st.name == "temperature"
        @test st.dims_xyzt == (true, true, true, true)
        @test st.nlayers == Nz
        @test length(st.depths) == Nz
        @test all(st.depths .> 0)               # positive down
        @test issorted(st.depths, rev = true)   # field order: deepest layer first
        @test st.unit == "°C"

        # 2D (surface) variable
        skin = TerrariumOutput(terrarium_model, :skin_temperature)
        @test skin.dims_xyzt == (true, true, false, true)
        @test skin.nlayers == 1

        # kwargs override derived metadata
        renamed = TerrariumOutput(terrarium_model, :temperature, name = "soil_temperature", keepbits = 10)
        @test renamed.name == "soil_temperature"
        @test renamed.keepbits == 10

        # unknown and unsupported (Face-located) variables throw
        @test_throws ArgumentError TerrariumOutput(terrarium_model, :not_a_variable)
        @test_throws ArgumentError TerrariumOutput(terrarium_model, :hydraulic_conductivity)

        # collect-all constructor: tuple of variables, no unsupported ones,
        # inputs excluded by default
        variables = TerrariumOutput(terrarium_model)
        @test variables isa Tuple
        names = [var.name for var in variables]
        @test "temperature" in names
        @test "skin_temperature" in names
        @test "saturation_water_ice" in names
        @test "hydraulic_conductivity" ∉ names
        @test "air_temperature" ∉ names    # input variable

        with_inputs = TerrariumOutput(terrarium_model, inputs = true)
        @test "air_temperature" in [var.name for var in with_inputs]
    end

    @testset "netCDF output run" begin
        land = SpeedyWeather.LandModel(spectral_grid, terrarium_model; Δt = 300.0)

        # default nlayers_soil: tests that the soil layer assertion in
        # initialize!(::NetCDFOutput, ...) is skipped without built-in land output
        tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
        output = NetCDFOutput(spectral_grid, PrimitiveWet, path = tmp_output_path, interval = Minute(30))

        model = PrimitiveWetModel(
            spectral_grid;
            land,
            land_sea_mask,
            output,
            surface_heat_flux = SurfaceHeatFlux(spectral_grid, land = PrescribedLandHeatFlux()),
            surface_humidity_flux = SurfaceHumidityFlux(spectral_grid, land = PrescribedLandHumidityFlux()),
            time_stepping = Leapfrog(spectral_grid, Δt_at_T31 = Minute(15)),
        )
        add!(model, TerrariumOutput(terrarium_model)...)

        simulation = SpeedyWeather.initialize!(model)
        SpeedyWeather.run!(simulation, period = Hour(1), output = true)

        NCDataset(SpeedyWeather.get_output_path(simulation)) do ds
            # shared soil depth dimension with layer-centre depths as coordinates
            @test ds.dim["soil_depth"] == Nz
            @test Vector(ds["soil_depth"][:]) ≈ TerrariumOutput(terrarium_model, :temperature).depths
            @test ds["soil_depth"].attrib["units"] == "m"
            @test ds["soil_depth"].attrib["positive"] == "down"

            # 3D variable on (lon, lat, soil_depth, time), 2D on (lon, lat, time)
            @test NCDatasets.dimnames(ds["temperature"]) == ("lon", "lat", "soil_depth", "time")
            @test NCDatasets.dimnames(ds["skin_temperature"]) == ("lon", "lat", "time")
            @test ds["temperature"].attrib["units"] == "°C"

            # finite over land, missing/NaN over ocean (fill value is NaN)
            isfilled(t) = ismissing(t) || isnan(t)
            temperature = ds["temperature"][:, :, :, end]
            @test any(t -> !isfilled(t) && isfinite(t), temperature)
            @test any(isfilled, temperature)
            skin_temperature = ds["skin_temperature"][:, :, end]
            @test any(t -> !isfilled(t) && isfinite(t), skin_temperature)
            @test any(isfilled, skin_temperature)
        end
    end

    @testset "Zarr output run" begin
        land = SpeedyWeather.LandModel(spectral_grid, terrarium_model; Δt = 300.0)

        tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
        output = ZarrOutput(spectral_grid, PrimitiveWet, path = tmp_output_path, interval = Minute(30))

        model = PrimitiveWetModel(
            spectral_grid;
            land,
            land_sea_mask,
            output,
            surface_heat_flux = SurfaceHeatFlux(spectral_grid, land = PrescribedLandHeatFlux()),
            surface_humidity_flux = SurfaceHumidityFlux(spectral_grid, land = PrescribedLandHumidityFlux()),
            time_stepping = Leapfrog(spectral_grid, Δt_at_T31 = Minute(15)),
        )
        add!(model, TerrariumOutput(terrarium_model)...)

        simulation = SpeedyWeather.initialize!(model)
        SpeedyWeather.run!(simulation, period = Hour(1), output = true)

        g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))

        # shared soil depth dimension with layer-centre depths as coordinates
        @test length(g["soil_depth"]) == Nz
        @test Vector(g["soil_depth"][:]) ≈ TerrariumOutput(terrarium_model, :temperature).depths
        @test g["soil_depth"].attrs["units"] == "m"
        @test g["soil_depth"].attrs["positive"] == "down"

        # 3D variable on (lon, lat, soil_depth, time), 2D on (lon, lat, time);
        # Zarr stores _ARRAY_DIMENSIONS in row-major (reverse of Julia) order
        @test g["temperature"].attrs["_ARRAY_DIMENSIONS"] == ["time", "soil_depth", "lat", "lon"]
        @test g["skin_temperature"].attrs["_ARRAY_DIMENSIONS"] == ["time", "lat", "lon"]
        @test g["temperature"].attrs["units"] == "°C"

        # finite over land, NaN (fill value) over ocean
        temperature = g["temperature"][:, :, :, end]
        @test any(isfinite, temperature)
        @test any(isnan, temperature)
        skin_temperature = g["skin_temperature"][:, :, end]
        @test any(isfinite, skin_temperature)
        @test any(isnan, skin_temperature)
    end
end
