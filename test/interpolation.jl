@testset "Interpolate constant field" begin
    npoints = 100
    
    for Grid in (   FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    HEALPix4Grid)
        
        for NF in (Float32,Float64)
        
            A = zeros(Grid{NF},8)           # something small
            c = randn(NF)
            A.data .= c                     # constant data

            θs = 180*rand(npoints) .- 90    # some latitudes in [-90˚,90˚N]
            λs = 360*rand(npoints)          # some longitudes in [0˚,360˚E]
            As = SpeedyWeather.interpolate(θs,λs,A)

            for a in As
                @test a ≈ c
            end
        end
    end
end