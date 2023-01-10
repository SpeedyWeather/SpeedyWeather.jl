@testset "Orographies" begin
    initialize_speedy(PrimitiveEquation,orography=NoOrography)
    initialize_speedy(PrimitiveEquation,orography=EarthOrography)
    initialize_speedy(PrimitiveEquation,orography=ZonalRidge)
end