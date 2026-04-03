@testset "Variable system for parameterizations" begin

    struct MyParam end

    MyParam(::SpectralGrid) = MyParam()

    SpeedyWeather.variables(::MyParam) = (
        ParameterizationVariable(:myvar_grid2d, SpeedyWeather.Grid2D(), desc = "My variable", units = "1"),
        ParameterizationVariable(:myvar_grid3d, SpeedyWeather.Grid3D(), desc = "My variable", units = "1"),
        PrognosticVariable(:myvar_spectral2d, SpeedyWeather.Spectral2D(), desc = "My variable", units = "1"),
        PrognosticVariable(:myvar_spectral3d, SpeedyWeather.Spectral3D(), desc = "My variable", units = "1"),
        PrognosticVariable(:mylandvar_spectral3d, SpeedyWeather.Spectral3D(), desc = "My variable", units = "1", namespace = :land),
        PrognosticVariable(:myoceanvar_spectral3d, SpeedyWeather.Spectral3D(), desc = "My variable", units = "1", namespace = :ocean),
    )

    # init model
    spectral_grid = SpectralGrid()

    # pass on to the model constructor ()
    model = PrimitiveWetModel(
        spectral_grid,
        parameterizations = (:custom_parameterization,),     # the parameterizations list defining order
        custom_parameterization = MyParam()
    )                  # MyParam() to be held by the model
    simulation = initialize!(model)

    # check that the variables are there and have the right dimension
    vars, model = SpeedyWeather.unpack(simulation)

    @test haskey(vars.parameterizations, :myvar_grid2d)
    @test ndims(vars.parameterizations.myvar_grid2d) == 1
    @test haskey(vars.parameterizations, :myvar_grid3d)
    @test ndims(vars.parameterizations.myvar_grid3d) == 2
    @test haskey(vars.prognostic, :myvar_spectral2d)
    @test ndims(vars.prognostic.myvar_spectral2d) == 1
    @test haskey(vars.prognostic, :myvar_spectral3d)
    @test ndims(vars.prognostic.myvar_spectral3d) == 2
    @test haskey(vars.prognostic.land, :mylandvar_spectral3d)
    @test ndims(vars.prognostic.land.mylandvar_spectral3d) == 2
    @test haskey(vars.prognostic.ocean, :myoceanvar_spectral3d)
    @test ndims(vars.prognostic.ocean.myoceanvar_spectral3d) == 2
end
