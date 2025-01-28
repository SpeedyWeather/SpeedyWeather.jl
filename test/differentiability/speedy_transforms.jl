# Tests for SpeedyTransforms

@testset "Differentiability: Complete Transform Enzyme" begin
    # make a high level finite difference test of the whole transform
    # can't use Enzyme or ChainRule Test tools for tests for that
    for (i_grid, grid_type) in enumerate(grid_types)

        spectral_grid = SpectralGrid(Grid=grid_type, trunc=10, nlayers=1, dealiasing=grid_dealiasing[i_grid])
        S = SpectralTransform(spectral_grid, one_more_degree=true)
        dS = deepcopy(S)

        if fd_tests[i_grid]
            
            # forwards 
            grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            dgrid = zero(grid)
            specs = zeros(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
            
            # seed
            dspecs = zero(specs)
            fill!(dspecs, 1+1im)

            autodiff(Reverse, transform!, Const, Duplicated(specs, dspecs), Duplicated(grid, dgrid), Duplicated(S, dS))

            # new seed
            dspecs2 = zero(specs)
            fill!(dspecs2, 1+1im)

            # finite difference comparision, seeded with a one adjoint to get the direct gradient
            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform(x, S), dspecs2, grid)
            @test isapprox(dgrid, fd_vjp[1])

            ## now backwards, as the input for spec we use the output of the forward transform

            fill!(dspecs,0)
            grid = zeros(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            dgrid = similar(grid)
            fill!(dgrid, 1)

            autodiff(Reverse, transform!, Const, Duplicated(grid, dgrid), Duplicated(specs, dspecs), Duplicated(S, dS))

            # new seed 
            dgrid2 = similar(grid)
            fill!(dgrid2, 1)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform(x, S), dgrid2, specs)

            @test isapprox(dspecs, fd_vjp[1])
        end 

        # test that d S^{-1}(S(x)) / dx = dx/dx = 1 (starting in both domains) 
        # this only holds for exact transforms, like Gaussian grids 

        # start with grid (but with a truncated one)
        function transform_identity!(x_out::AbstractGridArray{T}, x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
            x_SH = zeros(LowerTriangularArray{Complex{T}}, S.lmax+1, S.mmax+1, S.nlayers)
            transform!(x_SH, x, S)
            transform!(x_out, x_SH, S)
            return nothing
        end 

        function transform_identity(x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
            x_copy = deepcopy(x)
            transform_identity!(x_copy, x, S) 
            return x_copy
        end
        
        grid = rand(S.Grid{spectral_grid.NF}, S.nlat_half, S.nlayers)
        spec = transform(grid, S)

        grid = transform(spec, S)
        grid_out = zero(grid)
        
        transform_identity!(grid_out, grid, S)
        @test isapprox(grid, grid_out)

        dgrid = similar(grid)
        fill!(dgrid, 1)

        dgrid_out = zero(grid_out)

        autodiff(Reverse, transform_identity!, Const, Duplicated(grid_out, dgrid_out), Duplicated(grid, dgrid), Duplicated(S, dS))

        @test all(isapprox.(dgrid, 1)) 
        # TODO: previously this test was broken, with a version that directly mutates in-place. 
        # FD yields the same non-one values though. 
        # Not sure why. Do we use such things in our model? 
        #
        #function transform_identity!(x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
        #   x_SH = zeros(LowerTriangularArray{Complex{T}}, S.lmax+1, S.mmax+1, S.nlayers)
        #   transform!(x_SH, x, S)
        #   transform!(x, x_SH, S)
        #   return nothing
        #end 
        # The FD comparision passes, but it takes a long time to compute, so it's commented out. 
        #dgrid2 = similar(grid)
        #fill!(dgrid2, 1)
        #fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_identity(x, S), dgrid2, grid)
        #@test isapprox(dgrid, fd_vjp[1], rtol=0.01)

        # now start with spectral space, exclude for other grid because of https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/626
        if fd_tests[i_grid]

            function transform_identity!(x::LowerTriangularArray{Complex{T}}, S::SpectralTransform{T}) where T
                x_grid = zeros(S.Grid{T}, S.nlat_half, S.nlayers)
                transform!(x_grid, x, S)
                transform!(x, x_grid, S)
                return nothing
            end 

            spec = transform(grid, S)
            spec_copy = deepcopy(spec)
            transform_identity!(spec, S)
            @test isapprox(spec, spec_copy)

            dspec = similar(spec)
            fill!(dspec, 1+im)

            autodiff(Reverse, transform_identity!, Const, Duplicated(spec, dspec), Duplicated(S, dS))

            @test all(all.([isapprox.(dspec[il,1,:], 1) for il in 1:S.lmax+1])) # m = 0 => Im = 0

            for i in eachmatrix(dspec)
                @test all(isapprox.(dspec[:,i][S.lmax+2:end], 1+im)) 
            end 
        end
    end 
end

# Tests for all other spectral gradient functions
# We test that gradients are non-zero and identical with their finite difference
# When easiliy possible we check with the analytical formula as well 
@testset "Differentiability: Spectral Gradients Enzyme" begin 
    for (i_grid, grid_type) in enumerate(grid_types)

        if fd_tests[i_grid]

            spectral_grid = SpectralGrid(Grid=grid_type, trunc=10, nlayers=1, dealiasing=grid_dealiasing[i_grid])
            S = SpectralTransform(spectral_grid, one_more_degree=true)
            dS = deepcopy(S)

            u_grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            v_grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)

            u = transform(u_grid, S)
            v = transform(v_grid, S)
            du = zero(u)
            dv = zero(v)

            cu = zero(u)
            dcu = zero(u)
            fill!(dcu, 1)

            # curl test
            autodiff(Reverse, curl!, Const, Duplicated(cu, dcu), Duplicated(u, du), Duplicated(v, dv), Duplicated(S, dS))
            
            # we know  the gradient of the divergence wrt v is easy: im*m
            # so with a seed of 1, we should get for dv: im*m * 1 = im*m 
            # See https://speedyweather.github.io/SpeedyWeather.jl/dev/spectral_transform/
            # let's check it
            # TO-DO: why the other sign? but it's the same for Finite Differences
            # It's because it's the adjoint (')? And this matters here for complex numbers (see e.g. FiniteDifferences.jl examples for j'vp)
            # To-Do: double check that
            for i=1:dv.n
                @test all(Array(dv[:,1])[i:dv.m-1,i] .≈ complex(0,-(i-1)))
            end  
            @test sum(du) != 0 # nonzero gradient

            # new seed
            dcu2 = zero(dcu)
            fill!(dcu2, 1)

            # finite difference comparision, seeded with a one adjoint to get the direct gradient
            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> curl(x[1],x[2], S), dcu2, (u, v))
            @test isapprox(du, fd_vjp[1][1], rtol=1e-6)
            @test isapprox(dv, fd_vjp[1][2], rtol=1e-6)

            # div test
            u = transform(u_grid, S)
            v = transform(v_grid, S) 
            du = zero(u)
            dv = zero(v)
            div = zero(u)
            ddiv = zero(u)
            fill!(ddiv, 1 + 1im)

            autodiff(Reverse, divergence!, Const, Duplicated(div, ddiv), Duplicated(u, du), Duplicated(v, dv), Duplicated(S, dS))

            # we know the gradient of the divergence wrt u is easy: im*m 
            # See https://speedyweather.github.io/SpeedyWeather.jl/dev/spectral_transform/
            # let's check it 
            # To-Do: why the minus sign? 
            # It's because it's the adjoint (')? And this matters here for complex numbers
            # To-Do: double check that
            for i=1:du.n
                @test all(Array(du[:,1])[i:du.m-1,i] .≈ complex(i-1,-(i-1)))
            end 

            ddiv2 = zero(ddiv)
            fill!(ddiv2, 1 + 1im)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> divergence(x[1],x[2], S), ddiv2, (u, v))
            @test isapprox(du, fd_vjp[1][1])
            @test isapprox(dv, fd_vjp[1][2])
            @test sum(du) != 0 # nonzero gradient
            @test sum(dv) != 0 # nonzero gradient

            # UV_from_vor! 

            u = zero(u)
            du = fill!(du, 1+1im)

            v = zero(v)
            dv = fill!(dv, 1+1im)

            vor_grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            vor = transform(vor_grid, S)
            dvor = zero(vor)

            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.UV_from_vor!, Const, Duplicated(u, du), Duplicated(v, dv), Duplicated(vor, dvor), Duplicated(S, dS))

            function uvfvor(vor, S)
                u = zero(vor)
                v = zero(vor)
                SpeedyWeather.SpeedyTransforms.UV_from_vor!(u, v, vor, S)
                return cat(u, v, dims=2)
            end 
            
            uv_input = cat(u, v, dims=2)
            duv_input = zero(uv_input)
            duv_input = fill!(duv_input, 1+im)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> uvfvor(x, S), duv_input, vor)
            @test isapprox(dvor, fd_vjp[1])
            @test sum(dvor) != 0 # nonzero gradient

            # UV_from_vordiv! 
            u = zero(u)
            du = fill!(du, 1+1im)

            v = zero(v)
            dv = fill!(dv, 1+1im)

            vor_grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            vor = transform(vor_grid, S)
            dvor = zero(vor)

            div_grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            div = transform(vor_grid, S)
            ddiv = zero(vor)

            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.UV_from_vordiv!, Const, Duplicated(u, du), Duplicated(v, dv), Duplicated(vor, dvor), Duplicated(div, ddiv), Duplicated(S, dS))

            function uvfromvordiv(vor, div, S)
                u = zero(vor)
                v = zero(vor)
                SpeedyWeather.SpeedyTransforms.UV_from_vordiv!(u, v, vor, div, S)
                return cat(u, v, dims=2)
            end 
            
            uv_input = zero(uv_input)
            duv_input = fill!(duv_input, 1+im)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x-> uvfromvordiv(x[1], x[2], S), duv_input, (vor, div))
            @test isapprox(dvor, fd_vjp[1][1][:,1]) 
            @test isapprox(ddiv, fd_vjp[1][2][:,1])
            @test sum(dvor) != 0 # nonzero gradient
            @test sum(ddiv) != 0 # nonzero gradient

            # ∇²
            dvor = zero(vor)
            res_∇ = zero(vor)
            dres_∇ = zero(res_∇)
            fill!(dres_∇, 1+im)

            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.∇²!, Const, Duplicated(res_∇, dres_∇), Duplicated(vor, dvor), Duplicated(S, dS))

            dres_∇2 = zero(res_∇)
            fill!(dres_∇2, 1+im)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x-> ∇²(x, S), dres_∇2, vor)
            @test sum(dvor) != 0 # non-zero 
            @test isapprox(dvor, fd_vjp[1]) # and identical with FD

            # test with the eigenvalues saved in S, result should just be seed * eigenvalues
            for i=1:(vor.m-1)
                @test all(isapprox.(Array(dvor[:,1])[i,1:i], S.eigenvalues[i] * (1+im)))
            end 

            # ∇
            zonal_gradient = zero(vor)
            dzonal_gradient = zero(vor)
            fill!(dzonal_gradient, 1+im)

            merid_gradient = zero(vor)
            dmerid_gradient = zero(vor)
            fill!(dmerid_gradient, 1+im)

            dvor = zero(vor)
            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.∇!, Const, Duplicated(zonal_gradient, dzonal_gradient), Duplicated(merid_gradient, dmerid_gradient), Duplicated(vor, dvor), Duplicated(S, dS))

            dmerid_gradient2 = zero(dmerid_gradient)
            fill!(dmerid_gradient2, 1+im)

            dzonal_gradient2 = zero(dzonal_gradient)
            fill!(dzonal_gradient2, 1+im)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x-> ∇(x, S), (dmerid_gradient2, dzonal_gradient2), vor)
            @test sum(dvor) != # nonzero 
            @test isapprox(dvor, fd_vjp[1]) # and identical with FD
        end 
    end 
end 
 
