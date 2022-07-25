@testset "Model hierarchy" begin
    _,_,m_barotrop = initialize_speedy(model=:barotropic)
    _,_,m_shalloww = initialize_speedy(model=:shallowwater)
    _,_,m_primitive = initialize_speedy(model=:primitive)

    @test m_barotrop < m_shalloww
    @test m_shalloww < m_primitive
    
    @test (m_barotrop < m_barotrop) == false
    @test (m_shalloww < m_shalloww) == false
    @test (m_primitive < m_primitive) == false

    @test m_barotrop == m_barotrop
    @test m_shalloww == m_shalloww
    @test m_primitive == m_primitive

    # @test typeof(m_barotrop) <: SpeedyWeather.BarotropicModel
    # @test typeof(m_shalloww) <: SpeedyWeather.ShallowWaterModel
    # @test typeof(m_primitive) <: SpeedyWeather.PrimitiveEquationModel

    # @test SpeedyWeather.BarotropicModel < SpeedyWeather.ShallowWaterModel
    # @test SpeedyWeather.ShallowWaterModel < SpeedyWeather.PrimitiveEquationModel
    # @test SpeedyWeather.BarotropicModel < SpeedyWeather.PrimitiveEquationModel
end