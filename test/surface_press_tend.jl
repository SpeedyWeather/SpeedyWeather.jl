@testset "Surface pressure tendency no errors" begin
    for NF in (Float32,Float64)
        p,d,m = initialize_speedy(NF,model=PrimitiveEquation)
        SpeedyWeather.surface_pressure_tendency!(p,d,1,m)
    end
end