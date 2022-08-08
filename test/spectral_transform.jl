@testset "Transform: Triangular truncation" begin
    Ts = (31,42,85,170,341,682)     # spectral resolutions
    nlons = (96,128,256,512,1024)   # number of longitudes
    nlats = (48,64,128,256,512)     # number of latitudes
    for (T,nlon,nlat) in zip(Ts,nlons,nlats)
        
        tri_trunc = SpeedyWeather.triangular_truncation(trunc=T)    # trunc -> nlon, nlat
        @test (nlon,nlat) == (tri_trunc.nlon,tri_trunc.nlat)

        tri_trunc = SpeedyWeather.triangular_truncation(;nlon)      # nlon -> trunc
        @test T == tri_trunc.trunc
        
        tri_trunc = SpeedyWeather.triangular_truncation(;nlat)      # nlat -> trunc
        @test T == tri_trunc.trunc
    end
end

# for the following testsets test some spectral truncations
# but not too large ones as they take so long
spectral_resolutions = (31,42,85,170)

@testset "Transform: l=0,m=0 is constant > 0" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)

            p,d,m = initialize_speedy(NF;trunc)
            S = m.spectral_transform

            alms = copy(p.layers[1].leapfrog[1].vor)
            fill!(alms,0)
            alms[1,1] = 1

            map = gridded(alms,S)
        
            for i in eachindex(map)
                @test map[i] == map[1] > zero(NF)
            end
        end
    end
end

@testset "Transform: Recompute, precompute identical results" begin
    for trunc in spectral_resolutions
        for NF in (Float32,Float64)
            p1,d1,m1 = initialize_speedy(NF;trunc,recompute_legendre=false)
            p2,d2,m2 = initialize_speedy(NF;trunc,recompute_legendre=true)

            alms = copy(p1.layers[1].leapfrog[1].vor)
            alms .= randn(Complex{NF},size(alms)...)

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
            P = Parameters(;NF,trunc)
            S = SpectralTransform(P)

            lmax = 3
            for l in 1:lmax
                for m in 1:l
                    alms = zeros(LowerTriangularMatrix{Complex{NF}},S.lmax+2,S.mmax+1)
                    alms[l,m] = 1

                    map = gridded(alms,S)
                    alms2 = spectral(map,S)

                    for i in 1:S.lmax
                        for j in 1:i
                            @test alms[i,j] ≈ alms2[i,j] atol=100*eps(NF)
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
            P = Parameters(;NF,trunc)
            G = GeoSpectral(P)
            B = Boundaries(P)
            S = G.spectral_transform

            geopot_surf_spectral = B.geopot_surf
            geopot_surf_grid = gridded(geopot_surf_spectral,S)
            geopot_surf_spectral2 = spectral(geopot_surf_grid,S,one_more_l=true)
            SpeedyWeather.spectral_truncation!(geopot_surf_spectral2,trunc)
            geopot_surf_grid2 = gridded(geopot_surf_spectral2,S)

            for i in eachindex(geopot_surf_spectral)
                @test geopot_surf_spectral[i] ≈ geopot_surf_spectral2[i] rtol=30*sqrt(eps(NF))
            end
            for i in eachindex(geopot_surf_grid)
                @test geopot_surf_grid[i] ≈ geopot_surf_grid2[i] rtol=30*sqrt(eps(NF))
            end
        end
    end
end

@testset "Transform: with one more l" begin
    @testset for NF in (Float32,Float64)
        p,d,m = initialize_speedy(  NF,
                                    model=:shallowwater,
                                    initial_conditions=:rest,
                                    layer_thickness=0)

        # make sure vorticity and divergence are 0
        fill!(p.vor,0)
        fill!(p.div,0)

        # make sure vorticity and divergence are 0
        fill!(d.tendencies.vor_tend,0)                  
        fill!(d.tendencies.div_tend,0)

        # create initial conditions
        vor0 = zero(p.vor[:,:,1,1])
        div0 = zero(p.div[:,:,1,1])
        vor0[:,:] .= randn(Complex{NF},size(vor0)...)
        div0[:,:] .= randn(Complex{NF},size(div0)...)

        # remove non-zero entries in upper triangle
        SpeedyWeather.spectral_truncation!(vor0)
        SpeedyWeather.spectral_truncation!(div0)

        p.vor[:,:,1,1] .= vor0
        p.div[:,:,1,1] .= div0

        # get corresponding irrotational u_grid, v_grid (incl *coslat scaling)
        SpeedyWeather.gridded!(d,p,m)   

        # check we've actually created non-zero U,V
        @test all(d.grid_variables.U_grid .!= 0)
        @test all(d.grid_variables.V_grid .!= 0)

        lmax,mmax = size(vor0)      # 1-based maximum l,m
        
        @testset for lmax in (lmax,lmax+1)

            U = zeros(Complex{NF},lmax,mmax,1)
            U_grid = zero(d.grid_variables.U_grid)
            U_grid2 = zero(d.grid_variables.U_grid)

            # transform back and forth first as the original U_grid
            # is not exactly representable for lmax=lmax,lmax+1 (spectral truncation)
            S = m.geospectral.spectral_transform
            SpeedyWeather.spectral!(U,d.grid_variables.U_grid,S)
            SpeedyWeather.gridded!(U_grid,U,S)
            SpeedyWeather.spectral!(U,U_grid,S)
            SpeedyWeather.gridded!(U_grid2,U,S)

            testall = .≈(U_grid,U_grid2,rtol=sqrt(eps(NF)))
            @test sum(testall) > 0.99*length(testall)

            for i in eachindex(U_grid, U_grid2)
                @test U_grid[i] ≈ U_grid2[i] rtol=5*sqrt(sqrt(eps(NF)))
            end
        end
    end
end