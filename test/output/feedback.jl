@testset "Progress and parameter without output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_feedback_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid(nlayers=1)
    model = BarotropicModel(spectral_grid)
    simulation = initialize!(model)
    
    add!(model, ProgressTxt(path=tmp_output_path, write_only_with_output=true))
    add!(model, ParametersTxt(path=tmp_output_path, write_only_with_output=true))
    
    run!(simulation, period=Day(1))

    # test that files are not created because output=false
    @test ~isfile(joinpath(tmp_output_path, "parameters.txt"))
    @test ~isfile(joinpath(tmp_output_path, "progress.txt"))

    add!(model, ProgressTxt(path=tmp_output_path, write_only_with_output=false))
    add!(model, ParametersTxt(path=tmp_output_path, write_only_with_output=false))

    run!(simulation, period=Day(1))

    # test that files are created even if output=false
    @test isfile(joinpath(tmp_output_path, "parameters.txt"))
    @test isfile(joinpath(tmp_output_path, "progress.txt"))
end 