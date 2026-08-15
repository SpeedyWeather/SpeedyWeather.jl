# Tests for SpeedyTransforms
grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms
grid_dealiasing = [2.0, 3.0]
fd_tests = [true, true]
@testset "Differentiability: Complete Transform Enzyme" begin
    # make a high level finite difference test of the whole transform
    # can't use Enzyme or ChainRule Test tools for tests for that
    for (i_grid, grid_type) in enumerate(grid_types)

        spectral_grid = SpectralGrid(Grid = grid_type, truncation = 11, nlayers = 1, dealiasing = grid_dealiasing[i_grid])
        S = SpectralTransform(spectral_grid)
        dS = deepcopy(S)

        (; NF, grid, spectrum, nlayers) = spectral_grid

        if fd_tests[i_grid]

            # forwards
            field = rand(NF, grid, nlayers)
            dfield = zero(field)
            specs = zeros(complex(NF), spectrum, nlayers)

            # seed
            dspecs = zero(specs)
            fill!(dspecs, 1 + 1im)

            autodiff(Reverse, transform!, Const, Duplicated(specs, dspecs), Duplicated(field, dfield), Duplicated(S, dS))

            # new seed
            dspecs2 = zero(specs)
            fill!(dspecs2, 1 + 1im)

            # finite difference comparision, seeded with a one adjoint to get the direct gradient
            fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> transform(x, S), dspecs2, field)
            @test isapprox(dfield, fd_vjp[1])

            ## now backwards, as the input for spec we use the output of the forward transform

            fill!(dspecs, 0)
            field = zeros(NF, grid, nlayers)
            dfield = similar(field)
            fill!(dfield, 1)

            autodiff(Reverse, transform!, Const, Duplicated(field, dfield), Duplicated(specs, dspecs), Duplicated(S, dS))

            # new seed
            dfield2 = similar(field)
            fill!(dfield2, 1)

            fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> transform(x, S), dfield2, specs)

            @test isapprox(dspecs, fd_vjp[1])
        end

        # test that d S^{-1}(S(x)) / dx = dx/dx = 1 (starting in both domains)
        # this only holds for exact transforms, like Gaussian grids

        # start with field (but with a truncated one)
        # NOTE: size the intermediate to the LAYER COUNT OF `x`, not to `S.nlayers`. The latter is the
        # transform's scratch capacity (`max(maximum(transform_batch), 4*nlayers + 1)`, see the
        # `SpectralGrid` constructor in dynamics/spectral_grid.jl) and is deliberately larger than any
        # single call needs — using it here allocates a mismatched array and `transform!` throws.
        function transform_identity!(x_out::AbstractField, x::AbstractField, S::SpectralTransform{NF}) where {NF}
            x_SH = zeros(complex(NF), S.spectrum, size(x, 2))
            transform!(x_SH, x, S)
            transform!(x_out, x_SH, S)
            return nothing
        end

        function transform_identity(x::AbstractField, S::SpectralTransform)
            x_copy = deepcopy(x)
            transform_identity!(x_copy, x, S)
            return x_copy
        end

        field = rand(NF, grid, nlayers)
        spec = transform(field, S)

        field = transform(spec, S)
        field_out = zero(field)

        transform_identity!(field_out, field, S)
        @test isapprox(field, field_out)

        dfield = similar(field)
        fill!(dfield, 1)

        dfield_out = zero(field_out)

        autodiff(Reverse, transform_identity!, Const, Duplicated(field_out, dfield_out), Duplicated(field, dfield), Duplicated(S, dS))

        @test all(isapprox.(dfield, 1))
        # TODO: previously this test was broken, with a version that directly mutates in-place.
        # FD yields the same non-one values though.
        # Not sure why. Do we use such things in our model?
        #
        #function transform_identity!(x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
        #   x_SH = zeros(LowerTriangularArray{Complex{T}}, S.spectrum, S.nlayers)
        #   transform!(x_SH, x, S)
        #   transform!(x, x_SH, S)
        #   return nothing
        #end
        # The FD comparision passes, but it takes a long time to compute, so it's commented out.
        #dgrid2 = similar(grid)
        #fill!(dgrid2, 1)
        #fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_identity(x, S), dgrid2, grid)
        #@test isapprox(dfield, fd_vjp[1], rtol=0.01)

        # now start with spectral space, exclude for other grid because of https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/626
        if fd_tests[i_grid]

            function transform_identity!(x::LowerTriangularArray, S::SpectralTransform{NF}) where {NF}
                x_grid = zeros(NF, S.grid, size(x, 2))      # layer count of `x`, not S.nlayers (see above)
                transform!(x_grid, x, S)
                transform!(x, x_grid, S)
                return nothing
            end

            spec = transform(field, S)
            spec_copy = deepcopy(spec)
            transform_identity!(spec, S)
            @test isapprox(spec, spec_copy)

            dspec = similar(spec)
            fill!(dspec, 1 + im)

            autodiff(Reverse, transform_identity!, Const, Duplicated(spec, dspec), Duplicated(S, dS))

            @test all(all.([isapprox.(dspec[il, 1, :], 1) for il in 1:S.spectrum.lmax])) # m = 0 => Im = 0

            for i in eachmatrix(dspec)
                @test all(isapprox.(dspec[:, i][(S.spectrum.lmax + 1):end], 1 + im))
            end
        end
    end
end

# Tests for all other spectral gradient functions
# We test that gradients are non-zero and identical with their finite difference
# When easiliy possible we check with the analytical formula as well
# Real inner product over complex spectral arrays (treats them as real vectors), used for the
# adjoint identity below.
_realdot(a, b) = sum(real.(vec(Array(a))) .* real.(vec(Array(b)))) +
    sum(imag.(vec(Array(a))) .* imag.(vec(Array(b))))

# Tolerance for the adjoint identity below. These run at the `SpectralGrid`'s number format, which
# is Float32 by default, so a fixed 1e-10 is unreachable (Float32 eps is ~1.2e-7 and the dot
# products accumulate over all coefficients). `sqrt(eps)` is ~3.4e-4 in Float32 — still ~40x tighter
# than the ~1.4e-2 discrepancy a corrupted transform produced, so it retains its diagnostic power.
_adjoint_rtol(::Type{NF}) where {NF} = sqrt(eps(real(NF)))

# Every operator tested below (curl, divergence, ∇, ∇², UV_from_vor(div)) is LINEAR in its input, so
# applying the primal to a direction `v` already gives the Jacobian-vector product L·v. The adjoint
# identity  <L v, w> == <v, Lᵀ w>  then pins the reverse-mode result exactly, with no finite
# difference step size and no dependence on a complex-derivative convention.

@testset "Differentiability: Spectral Gradients Enzyme" begin
    for (i_grid, grid_type) in enumerate(grid_types)

        if fd_tests[i_grid]

            spectral_grid = SpectralGrid(Grid = grid_type, truncation = 11, nlayers = 1, dealiasing = grid_dealiasing[i_grid])

            # Reverse-mode autodiff MUTATES the primal transform: it accumulates into
            # `S.gradients.grad_y_vordiv1/2`, even though `SpectralTransform` is declared
            # `EnzymeRules.inactive_type` and even when it is passed as `Const`. Forward mode and
            # plain (non-AD) calls leave it alone. Reusing one transform across several `autodiff`
            # calls therefore silently computes every call after the first with a corrupted
            # operator — which is what made the checks below disagree with the analytic results.
            # So: a fresh transform per `autodiff` call, plus `S_ref`, which is never handed to
            # `autodiff` and is used for all primal/JVP references.
            fresh_transform() = SpectralTransform(spectral_grid, one_more_degree = true)
            S_ref = fresh_transform()
            S = fresh_transform()
            dS = deepcopy(S)

            (; NF, grid, spectrum, nlayers) = spectral_grid
            u_grid = rand(NF, grid, nlayers)
            v_grid = rand(NF, grid, nlayers)

            u = transform(u_grid, S_ref)
            v = transform(v_grid, S_ref)
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
            for i in 1:dv.spectrum.mmax
                @test all(Array(dv[:, 1])[i:(dv.spectrum.lmax - 1), i] .≈ complex(0, -(i - 1)))
            end
            @test sum(du) != 0 # nonzero gradient

            # adjoint identity for curl: <curl(vu, vv), w> == <vu, du> + <vv, dv> with w the seed
            # used above (`dcu` before it was consumed, i.e. all ones)
            let vu = rand(complex(NF), spectrum, nlayers), vv = rand(complex(NF), spectrum, nlayers)
                w = zero(cu); fill!(w, 1)
                Lv = zero(cu); curl!(Lv, vu, vv, S_ref)      # linear ⇒ the primal IS the JVP
                @test isapprox(_realdot(Lv, w), _realdot(vu, du) + _realdot(vv, dv); rtol = _adjoint_rtol(NF))
            end

            # div test
            u = transform(u_grid, S_ref)
            v = transform(v_grid, S_ref)
            du = zero(u)
            dv = zero(v)
            div = zero(u)
            ddiv = zero(u)
            fill!(ddiv, 1 + 1im)

            S = fresh_transform(); dS = deepcopy(S)
            autodiff(Reverse, divergence!, Const, Duplicated(div, ddiv), Duplicated(u, du), Duplicated(v, dv), Duplicated(S, dS))

            # we know the gradient of the divergence wrt u is easy: im*m
            # See https://speedyweather.github.io/SpeedyWeather.jl/dev/spectral_transform/
            # let's check it
            # To-Do: why the minus sign?
            # It's because it's the adjoint (')? And this matters here for complex numbers
            # To-Do: double check that
            for i in 1:du.spectrum.mmax
                @test all(Array(du[:, 1])[i:(du.spectrum.lmax - 1), i] .≈ complex(i - 1, -(i - 1)))
            end

            # adjoint identity for divergence, with the same seed (1 + 1im) as used above
            let vu = rand(complex(NF), spectrum, nlayers), vv = rand(complex(NF), spectrum, nlayers)
                w = zero(div); fill!(w, 1 + 1im)
                Lv = zero(div); divergence!(Lv, vu, vv, S_ref)
                @test isapprox(_realdot(Lv, w), _realdot(vu, du) + _realdot(vv, dv); rtol = _adjoint_rtol(NF))
            end
            @test sum(du) != 0 # nonzero gradient
            @test sum(dv) != 0 # nonzero gradient

            # UV_from_vor!

            u = zero(u)
            du = fill!(du, 1 + 1im)

            v = zero(v)
            dv = fill!(dv, 1 + 1im)

            vor_grid = rand(NF, spectral_grid.grid, spectral_grid.nlayers)
            vor = transform(vor_grid, S_ref)
            dvor = zero(vor)

            S = fresh_transform(); dS = deepcopy(S)
            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.UV_from_vor!, Const, Duplicated(u, du), Duplicated(v, dv), Duplicated(vor, dvor), Duplicated(S, dS))

            # adjoint identity for UV_from_vor! (one input, two outputs), seeds as used above
            let vvor = rand(complex(NF), spectrum, nlayers)
                wu = zero(u); fill!(wu, 1 + 1im)
                wv = zero(v); fill!(wv, 1 + 1im)
                Lu = zero(u); Lv = zero(v)
                SpeedyWeather.SpeedyTransforms.UV_from_vor!(Lu, Lv, vvor, S_ref)
                @test isapprox(_realdot(Lu, wu) + _realdot(Lv, wv), _realdot(vvor, dvor); rtol = _adjoint_rtol(NF))
            end
            @test sum(dvor) != 0 # nonzero gradient

            # UV_from_vordiv!
            u = zero(u)
            du = fill!(du, 1 + 1im)

            v = zero(v)
            dv = fill!(dv, 1 + 1im)

            vor_grid = rand(NF, grid, spectral_grid.nlayers)
            vor = transform(vor_grid, S_ref)
            dvor = zero(vor)

            div_grid = rand(NF, grid, spectral_grid.nlayers)
            div = transform(div_grid, S_ref)
            ddiv = zero(vor)

            S = fresh_transform(); dS = deepcopy(S)
            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.UV_from_vordiv!, Const, Duplicated(u, du), Duplicated(v, dv), Duplicated(vor, dvor), Duplicated(div, ddiv), Duplicated(S, dS))

            # adjoint identity for UV_from_vordiv! (two inputs, two outputs)
            let vvor = rand(complex(NF), spectrum, nlayers), vdiv = rand(complex(NF), spectrum, nlayers)
                wu = zero(u); fill!(wu, 1 + 1im)
                wv = zero(v); fill!(wv, 1 + 1im)
                Lu = zero(u); Lv = zero(v)
                SpeedyWeather.SpeedyTransforms.UV_from_vordiv!(Lu, Lv, vvor, vdiv, S_ref)
                @test isapprox(
                    _realdot(Lu, wu) + _realdot(Lv, wv),
                    _realdot(vvor, dvor) + _realdot(vdiv, ddiv); rtol = _adjoint_rtol(NF),
                )
            end
            @test sum(dvor) != 0 # nonzero gradient
            @test sum(ddiv) != 0 # nonzero gradient

            # ∇²
            dvor = zero(vor)
            res_∇ = zero(vor)
            dres_∇ = zero(res_∇)
            fill!(dres_∇, 1 + im)

            S = fresh_transform(); dS = deepcopy(S)
            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.∇²!, Const, Duplicated(res_∇, dres_∇), Duplicated(vor, dvor), Duplicated(S, dS))

            # adjoint identity for ∇²
            let vvor = rand(complex(NF), spectrum, nlayers)
                w = zero(res_∇); fill!(w, 1 + 1im)
                Lv = zero(res_∇)
                SpeedyWeather.SpeedyTransforms.∇²!(Lv, vvor, S_ref)
                @test isapprox(_realdot(Lv, w), _realdot(vvor, dvor); rtol = _adjoint_rtol(NF))
            end
            @test sum(dvor) != 0 # non-zero

            # test with the eigenvalues saved in S, result should just be seed * eigenvalues
            # (eigenvalues live directly on the transform, not in `S.gradients` — see the `Gradients` struct)
            for i in 1:(vor.spectrum.mmax)
                @test all(isapprox.(Array(dvor[:, 1])[i, 1:i], S_ref.eigenvalues[i] * (1 + im)))
            end

            # ∇
            zonal_gradient = zero(vor)
            dzonal_gradient = zero(vor)
            fill!(dzonal_gradient, 1 + im)

            merid_gradient = zero(vor)
            dmerid_gradient = zero(vor)
            fill!(dmerid_gradient, 1 + im)

            dvor = zero(vor)
            S = fresh_transform(); dS = deepcopy(S)
            autodiff(Reverse, SpeedyWeather.SpeedyTransforms.∇!, Const, Duplicated(zonal_gradient, dzonal_gradient), Duplicated(merid_gradient, dmerid_gradient), Duplicated(vor, dvor), Duplicated(S, dS))

            # adjoint identity for ∇ (one input, two outputs: zonal and meridional gradient).
            # NOTE: the previous version read `@test sum(dvor) != # nonzero` followed by a second
            # `@test` on the next line — the `!=` swallowed it, so neither assertion tested what it
            # looked like it did.
            let vvor = rand(complex(NF), spectrum, nlayers)
                wz = zero(zonal_gradient); fill!(wz, 1 + 1im)
                wm = zero(merid_gradient); fill!(wm, 1 + 1im)
                Lz = zero(zonal_gradient); Lm = zero(merid_gradient)
                SpeedyWeather.SpeedyTransforms.∇!(Lz, Lm, vvor, S_ref)
                @test isapprox(_realdot(Lz, wz) + _realdot(Lm, wm), _realdot(vvor, dvor); rtol = _adjoint_rtol(NF))
            end
            @test sum(dvor) != 0 # nonzero
        end
    end
end
