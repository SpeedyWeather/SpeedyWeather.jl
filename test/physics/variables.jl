using Adapt

@testset "Test the variable system for the parameterization" begin
    
    struct MyParam end 
    
    MyParam(::SpectralGrid) = MyParam()

    Adapt.@adapt_structure MyParam

    SpeedyWeather.variables(::MyParam) = (DiagnosticVariable(name=:myvar_grid2d, dims=Grid2D(), desc="My variable", units="1"),
    
    DiagnosticVariable(name=:myvar_grid3d, dims=Grid3D(), desc="My variable", units="1"),
    PrognosticVariable(name=:myvar_spectral2d, dims=Spectral2D(), desc="My variable", units="1"),
    PrognosticVariable(name=:myvar_spectral3d, dims=Spectral3D(), desc="My variable", units="1"), 
    PrognosticVariable(name=:mylandvar_spectral3d, dims=Spectral3D(), desc="My variable", units="1", namespace=:land),
    DiagnosticVariable(name=:myoceanvar_spectral3d, dims=Spectral3D(), desc="My variable", units="1", namespace=:ocean))

    # init model 
    spectral_grid = SpectralGrid()

    # pass on to the model constructor ()
    model = PrimitiveWetModel(spectral_grid,
        parameterizations=(:custom_parameterization, ),     # the parameterizations list defining order
        custom_parameterization=MyParam())                  # MyParam() to be held by the model
    simulation = initialize!(model)

    # check that the variables are there and have the right dimension
    progn, diagn, model = SpeedyWeather.unpack(simulation)
    
    @test haskey(diagn.physics, :myvar_grid2d)
    @test ndims(diagn.physics.myvar_grid2d) == 1
    @test haskey(diagn.physics, :myvar_grid3d)
    @test ndims(diagn.physics.myvar_grid3d) == 2
    @test haskey(progn.physics, :myvar_spectral2d)
    @test ndims(progn.physics.myvar_spectral2d) == 1
    @test haskey(progn.physics, :myvar_spectral3d)
    @test ndims(progn.physics.myvar_spectral3d) == 2
    @test haskey(progn.land, :mylandvar_spectral3d)
    @test ndims(progn.land.mylandvar_spectral3d) == 2
    @test haskey(diagn.physics.ocean, :myoceanvar_spectral3d)
    @test ndims(diagn.physics.ocean.myoceanvar_spectral3d) == 2
end