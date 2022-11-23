import NetCDF: NetCDF

@testset "Download input data" begin
    input_data_path = normpath(joinpath(@__FILE__, "..", "..", "input_data"))
    rm(input_data_path, force=true, recursive=true)
    isdir(input_data_path) || SpeedyWeather.download_input_data(input_data_path);
    ncfile = NetCDF.open("../input_data/orography_F512.nc")
end
