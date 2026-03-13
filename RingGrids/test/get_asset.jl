using NCDatasets

# path to a known asset file in SpeedyWeatherAssets
const LSM_PATH = joinpath("data", "boundary_conditions", "land-sea_mask.nc")

@testset "get_asset: Download and load into Array (NCDatasets)" begin
    # load land-sea mask as Array with explicit variable name
    data = get_asset(LSM_PATH; name = "lsm", FileFormat = NCDataset)
    @test data isa Array
    @test ndims(data) == 2
    @test size(data, 1) > 0
    @test size(data, 2) > 0
    # land-sea mask values should be in [0, 1]
    @test all(0 .<= data .<= 1)
end

@testset "get_asset: Download and load into FullClenshawField (NCDatasets)" begin
    data = get_asset(LSM_PATH; name = "lsm", ArrayType = FullClenshawField, FileFormat = NCDataset)
    @test data isa FullClenshawField
    @test length(data) > 0
    @test all(0 .<= data .<= 1)
end

@testset "get_asset: Version keyword (NCDatasets)" begin
    # load with explicit version number
    data_v = get_asset(LSM_PATH; name = "lsm", FileFormat = NCDataset, version = v"1")
    @test data_v isa Array
    @test ndims(data_v) == 2

    # load from branch name
    data_branch = get_asset(LSM_PATH; name = "lsm", FileFormat = NCDataset, version = "main")
    @test data_branch isa Array
    @test size(data_branch) == size(data_v)
end

@testset "get_asset: Local file loading with NCDatasets (from_assets=false)" begin
    # first download the asset to get a local path
    data_remote = get_asset(LSM_PATH; name = "lsm", FileFormat = NCDataset)

    # find the local artifact path by going through the artifact system
    artifact_toml = joinpath(pkgdir(RingGrids), "Artifacts.toml")
    hash = Pkg.Artifacts.artifact_hash(basename(LSM_PATH), artifact_toml)
    local_path = joinpath(Pkg.Artifacts.artifact_path(hash), basename(LSM_PATH))

    # load the same file locally
    data_local = get_asset(local_path; from_assets = false, name = "lsm", FileFormat = NCDataset)
    @test data_local isa Array
    @test data_local == data_remote
end

@testset "get_asset: Nonexistent variable name (NCDatasets)" begin
    @test_throws ErrorException get_asset(
        LSM_PATH; name = "nonexistent_var", FileFormat = NCDataset
    )
end

@testset "get_asset: Interpolation to output_grid (NCDatasets)" begin
    # use orography (Float32) which is on a FullGaussianGrid
    oro_path = joinpath("data", "boundary_conditions", "orography.nc")
    output_grid = FullGaussianGrid(16)
    data = get_asset(
        oro_path; name = "orog", ArrayType = FullGaussianField,
        FileFormat = NCDataset, output_grid = output_grid
    )
    @test data isa FullGaussianField
    @test RingGrids.get_nlat_half(data.grid) == 16
    @test length(data) > 0
end
