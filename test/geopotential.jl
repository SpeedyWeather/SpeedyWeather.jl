@testset "Geopotential no errors" begin
    for NF in (Float32,Float64)
        p,d,m = initialize_speedy(NF,model=PrimitiveEquation)
        SpeedyWeather.geopotential!(d,p,1,m.boundaries,m.geometry)
    end
end