@testset "FullGaussianGrid: Test grid and spectral resolution match" begin
    Ts = (31,42,85,170,341,682)     # spectral resolutions
    nlons = (96,128,256,512,1024)   # number of longitudes
    nlats = (48,64,128,256,512)     # number of latitudes
    for (T,nlon,nlat) in zip(Ts,nlons,nlats)
        
        nlat_half = SpeedyWeather.get_resolution(FullGaussianGrid,T)
        @test (nlon,nlat) == (4nlat_half,2nlat_half)

        trunc = SpeedyWeather.get_truncation(FullGaussianGrid,nlat÷2)
        @test T == trunc
    end
end

# for the following testsets test some spectral truncations
# but not too large ones as they take so long
spectral_resolutions = (31,63,127)

@testset "Transform: l=0,m=0 is constant > 0" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,
                            HEALPixGrid)

                p,d,m = initialize_speedy(NF;trunc,Grid)
                S = m.spectral_transform

                alms = copy(p.layers[1].leapfrog[1].vor)
                fill!(alms,0)
                alms[1,1] = 1

                map = gridded(alms,S)
            
                for ij in SpeedyWeather.eachgridpoint(map)
                    @test map[ij] ≈ map[1] > zero(NF)
                end
            end
        end
    end
end

@testset "Transform: Recompute, precompute identical results" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            p1,d1,m1 = initialize_speedy(NF;trunc,recompute_legendre=false)
            p2,d2,m2 = initialize_speedy(NF;trunc,recompute_legendre=true)

            (;vor) = p1.layers[1].leapfrog[1]
            alms = randn(typeof(vor),size(vor)...)

            map1 = gridded(alms,m1.spectral_transform)
            map2 = gridded(alms,m2.spectral_transform)
        
            # is only approx as recompute_legendre may use a different precision
            @test map1 ≈ map2
        end
    end
end

@testset "Transform: Individual Legendre polynomials" begin
    @testset for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid)
                            # HEALPixGrid)
                P = Parameters(;NF,trunc,Grid)
                S = SpectralTransform(P)

                lmax = 3
                for l in 1:lmax
                    for m in 1:l
                        alms = zeros(LowerTriangularMatrix{Complex{NF}},S.lmax+2,S.mmax+1)
                        alms[l,m] = 1

                        map = gridded(alms,S)
                        alms2 = spectral(map,S)

                        for lm in SpeedyWeather.eachharmonic(alms,alms2)
                            @test alms[lm] ≈ alms2[lm] atol=100*eps(NF)
                        end
                    end
                end
            end
        end
    end
end

@testset "Transform: Geopotential" begin

    # Test for variable resolution
    @testset for trunc in spectral_resolutions[1:2]
        for NF in (Float64,Float32)
            for Grid in (   FullGaussianGrid,
                            FullClenshawGrid,
                            OctahedralGaussianGrid,
                            OctahedralClenshawGrid,)
                            # HEALPixGrid)
                P = Parameters(;NF,trunc,model=:shallowwater)
                S = SpectralTransform(P)
                B = Boundaries(P,S)

                oro_grid = B.orography
                oro_spec = spectral(oro_grid,S)
                oro_grid1 = gridded(oro_spec,S)
                oro_spec1 = spectral(oro_grid1,S)
                oro_grid2 = gridded(oro_spec1,S)
                oro_spec2 = spectral(oro_grid2,S)

                for lm in SpeedyWeather.eachharmonic(oro_spec1,oro_spec2)
                    @test oro_grid1[lm] ≈ oro_grid2[lm] rtol=200*sqrt(eps(NF))
                end
                for ij in eachindex(oro_grid1,oro_grid2)
                    @test oro_grid1[ij] ≈ oro_grid2[ij] rtol=200*sqrt(eps(NF))
                end
            end
        end
    end
end

# @testset "Transform: with one more l" begin
#     @testset for NF in (Float32,Float64)
#         p,d,m = initialize_speedy(  NF,
#                                     model=:shallowwater,
#                                     initial_conditions=:rest,
#                                     layer_thickness=0)

#         # make sure vorticity and divergence are 0
#         fill!(p.vor,0)
#         fill!(p.div,0)

#         # make sure vorticity and divergence are 0
#         fill!(d.tendencies.vor_tend,0)                  
#         fill!(d.tendencies.div_tend,0)

#         # create initial conditions
#         vor0 = zero(p.vor[:,:,1,1])
#         div0 = zero(p.div[:,:,1,1])
#         vor0[:,:] .= randn(Complex{NF},size(vor0)...)
#         div0[:,:] .= randn(Complex{NF},size(div0)...)

#         # remove non-zero entries in upper triangle
#         SpeedyWeather.spectral_truncation!(vor0)
#         SpeedyWeather.spectral_truncation!(div0)

#         p.vor[:,:,1,1] .= vor0
#         p.div[:,:,1,1] .= div0

#         # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
#         SpeedyWeather.gridded!(d,p,m)   

#         # check we've actually created non-zero U,V
#         @test all(d.grid_variables.U_grid .!= 0)
#         @test all(d.grid_variables.V_grid .!= 0)

#         lmax,mmax = size(vor0)      # 1-based maximum l,m
        
#         @testset for lmax in (lmax,lmax+1)

#             U = zeros(Complex{NF},lmax,mmax,1)
#             U_grid = zero(d.grid_variables.U_grid)
#             U_grid2 = zero(d.grid_variables.U_grid)

#             # transform back and forth first as the original U_grid
#             # is not exactly representable for lmax=lmax,lmax+1 (spectral truncation)
#             S = m.geospectral.spectral_transform
#             SpeedyWeather.spectral!(U,d.grid_variables.U_grid,S)
#             SpeedyWeather.gridded!(U_grid,U,S)
#             SpeedyWeather.spectral!(U,U_grid,S)
#             SpeedyWeather.gridded!(U_grid2,U,S)

#             testall = .≈(U_grid,U_grid2,rtol=sqrt(eps(NF)))
#             @test sum(testall) > 0.99*length(testall)

#             for i in eachindex(U_grid, U_grid2)
#                 @test U_grid[i] ≈ U_grid2[i] rtol=5*sqrt(sqrt(eps(NF)))
#             end
#         end
#     end
# end
