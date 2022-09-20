@testset "Initialize sigma levels manually" begin

    # automatic levels
    p,d,m = initialize_speedy(nlev=4)
    @test length(m.geometry.σ_levels_half) == 5
    @test length(m.geometry.σ_levels_full) == 4

    # manual levels
    p,d,m = initialize_speedy(σ_levels_half=[0,0.4,0.6,1])
    @test m.parameters.nlev == 3
    @test length(m.geometry.σ_levels_half) == 4
    @test length(m.geometry.σ_levels_full) == 3

    # specify both 
    p,d,m = initialize_speedy(σ_levels_half=[0,0.4,0.6,1],nlev=3)
    @test m.parameters.nlev == 3
    @test length(m.geometry.σ_levels_half) == 4
    @test length(m.geometry.σ_levels_full) == 3
      
end