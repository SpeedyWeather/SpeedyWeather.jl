import Random

@testset "Interpolate constant field" begin
    npoints = 100
    
    @testset for Grid in (   FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid)
        
        @testset for NF in (Float32, Float64)
        
            A = zeros(Grid{NF}, 8)           # something small
            c = randn(NF)
            A.data .= c                     # constant data

            θs = 180*rand(npoints) .- 90    # some latitudes in [-90˚, 90˚N]
            λs = 360*rand(npoints)          # some longitudes in [0˚, 360˚E]
            As = RingGrids.interpolate(θs, λs, A)

            for a in As
                @test a ≈ c
            end
        end
    end
end

@testset "Interpolate zonally-constant field" begin
    npoints = 1000
    
    @testset for Grid in (   FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid)
        
        @testset for NF in (Float32, Float64)
        
            A = zeros(Grid{NF}, 32)          # some resolution
            G = RingGrids.GridGeometry(A)
            lat1 = G.latd[2]                # latitude of first ring
            
            for (j, ring) in enumerate(RingGrids.eachring(A))
                θ = G.latd[j+1]     # G.latd also includes 90˚N hence +1
                for ij in ring
                    A[ij] = θ
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            θs = 2lat1*rand(npoints) .- lat1    # some latitudes in [-90˚, 90˚N]
            λs = 360*rand(npoints)              # some longitudes in [0˚, 360˚E]
            
            As = RingGrids.interpolate(θs, λs, A)

            for (a, θ) in zip(As, θs)
                @test a ≈ θ
            end
        end
    end
end

@testset "Interpolate meridionally-constant field" begin
    npoints = 1000
    
    @testset for Grid in (  FullGaussianGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,
                            HEALPixGrid,
                            OctaHEALPixGrid)
        
        @testset for NF in (Float32, Float64)
        
            A = zeros(Grid{NF}, 32)          # some resolution
            G = RingGrids.GridGeometry(A)
            lat1 = G.latd[2]                # latitude of first ring
            
            # TEST FROM 60˚ to 300˚ to not interpolate across 0/360˚E
            # where this test won't work because lon have a sharp jump
            # and aren't linear across the prime meridian
            # but that differently further down
            for (j, ring) in enumerate(RingGrids.eachring(A))
                for ij in ring
                    A[ij] = G.londs[ij]
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            θs = 2lat1*rand(npoints) .- lat1    # some latitudes in (-90˚, 90˚N)
            λs = 240*rand(npoints) .+ 60        # some longitudes in [60˚, 300˚E]
            
            As = RingGrids.interpolate(θs, λs, A)

            for (a, λ) in zip(As, λs)
                @test a ≈ λ
            end

            f(λ) = λ > 180 ? λ-360 : λ          # 0-360˚ to -180˚-180˚E

            # TEST FROM -120˚ to 120˚ to still test the indexing across the
            # prime meridian
            for (j, ring) in enumerate(RingGrids.eachring(A))
                for ij in ring
                    A[ij] = f(G.londs[ij])
                end
            end

            # don't interpolate above first or below last ring
            # northpole value isn't 90 as it's just the average of the
            # first ring, same for south pole
            θs = 2lat1*rand(npoints) .- lat1    # some latitudes in (-90˚, 90˚N)
            λs = 240*rand(npoints) .- 120       # some longitudes in [-120˚, 120˚E]
            
            As = RingGrids.interpolate(θs, λs, A)

            for (a, λ) in zip(As, λs)
                @test a ≈ λ
            end
        end
    end
end

@testset "Find latitude rings and weights" begin
    @testset for Grid in (   FullGaussianGrid,
                    OctahedralGaussianGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaHEALPixGrid)

        @testset for nlat_half in [4, 8, 16] 
            G = RingGrids.GridGeometry(Grid, nlat_half)
            latd = G.latd

            n = length(latd)-1
            Δs = rand(n)

            r = Random.randperm(n)
            θs = latd[1:end-1] .+ diff(latd).*Δs

            θs = θs[r]
            Δs = Δs[r]
            js, Δys = RingGrids.find_rings(θs, latd)

            for (i, (j, Δref, Δ)) in enumerate(zip(js, Δs, Δys))
                @test j == r[i]-1
                @test Δref ≈ Δ
            end
        end
    end
end

@testset "Interpolate between grids" begin
    @testset for NF in (Float32, Float64)
        @testset for Grid in (  FullGaussianGrid,
                                FullClenshawGrid,
                                OctahedralGaussianGrid,
                                OctahedralClenshawGrid,
                                HEALPixGrid,
                                OctaHEALPixGrid)

            # create some smooth gridded field
            trunc = 10
            alms = randn(LowerTriangularMatrix{Complex{NF}}, 5, 5)
            alms = spectral_truncation(alms, trunc+2, trunc+1)
            A = transform(alms; Grid)

            # interpolate to FullGaussianGrid and back and compare
            nlat_half = 32
            A_interpolated = RingGrids.interpolate(FullGaussianGrid, nlat_half, A)
            A2 = zero(A)
            RingGrids.interpolate!(A2, A_interpolated)

            # just check that it's not completely off
            for ij in RingGrids.eachgridpoint(A, A2)
                @test A[ij] ≈ A2[ij] rtol=5e-1 atol=5e-1
            end
        end
    end
end