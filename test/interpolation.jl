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

# @testset "Interpolate zonally-constant field" begin
#     npoints = 4
    
#     for Grid in (   FullGaussianGrid,
#                     OctahedralGaussianGrid,
#                     OctahedralClenshawGrid,
#                     HEALPixGrid,
#                     HEALPix4Grid)
        
#         for NF in (Float32,Float64)
        
#             A = zeros(Grid{NF},8)           # something small
#             G = SpeedyWeather.GridGeometry(A)
            
#             for (j,ring) in enumerate(SpeedyWeather.eachring(A))
#                 θ = G.latd[j+1]     # G.latd also includes 90˚N
#                 for ij in ring
#                     A[ij] = θ
#                 end
#             end

#             θs = 180*rand(npoints) .- 90    # some latitudes in [-90˚,90˚N]
#             λs = 360*rand(npoints)          # some longitudes in [0˚,360˚E]
            
#             println(θs)
#             As = SpeedyWeather.interpolate(θs,λs,A)
#             println(As)
#             println("---")


#             # for (a,θ) in (As,θs)
#             #     @test a ≈ θ
#             # end
#         end
#     end
# end

@testset "Find latitude rings and weights" begin
    for Grid in (   FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    HEALPix4Grid)

        for nlat_half in [4,8,16] 
            G = SpeedyWeather.GridGeometry(Grid,nlat_half)
            latd = G.latd

            n = length(latd)-1
            Δs = rand(n)

            r = Random.randperm(n)
            θs = latd[1:end-1] .+ diff(latd).*Δs

            θs = θs[r]
            Δs = Δs[r]
            js,Δys = SpeedyWeather.find_rings(θs,latd)

            for (i,(j,Δref,Δ)) in enumerate(zip(js,Δs,Δys))
                @test j == r[i]-1
                @test Δref ≈ Δ
            end
        end
    end
end