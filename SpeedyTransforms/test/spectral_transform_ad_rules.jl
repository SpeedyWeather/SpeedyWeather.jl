using EnzymeTestUtils, Enzyme
import EnzymeTestUtils: test_approx
import AbstractFFTs
using FiniteDifferences

grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms
grid_dealiasing = [2, 3]
fd_tests = [true, true]

# currently there's an issue with EnzymeTestUtils not being able to work with structs with undefined fields like FFT plans
# https://github.com/EnzymeAD/Enzyme.jl/issues/1992
# This is a very hacky workaround
function EnzymeTestUtils.test_approx(x::AbstractFFTs.Plan, y::AbstractFFTs.Plan, msg; kwargs...)
    EnzymeTestUtils.@test_msg "$msg: types must match" typeof(x) == typeof(y)
    names = fieldnames(typeof(x))[1:(end - 1)] # exclude pinv field (which is the last field)
    if isempty(names)
        EnzymeTestUtils.@test_msg msg x == y
    else
        for k in names
            if k isa Symbol && hasproperty(x, k)
                msg_new = "$msg: ::$(typeof(x)).$k"
            else
                msg_new = "$msg: getfield(::$(typeof(x)), $k)"
            end
            EnzymeTestUtils.test_approx(getfield(x, k), getfield(y, k), msg_new; kwargs...)
        end
    end
    return nothing
end

# `test_forward` requires a mutating function to RETURN the argument(s) it mutates — the opposite of
# `test_reverse`, which needs it to return nothing. `_fourier!` returns nothing, so wrap it; the
# grid -> Fourier direction writes both scratch halves, so both are returned to get both tangents checked.
fourier_analysis!(f_north, f_south, field, S) =
    (SpeedyTransforms._fourier!(f_north, f_south, field, S); (f_north, f_south))
fourier_synthesis!(field, f_north, f_south, S) =
    (SpeedyTransforms._fourier!(field, f_north, f_south, S); field)

@testset "SpeedyTransforms: AD Rules" begin
    @testset "_fourier! Enzyme rules" begin
        @testset "EnzymeTestUtils reverse rule test" begin
            for (i_grid, Grid) in enumerate(grid_types)

                # these tests don't pass for reduced grids
                # this is likely due to FiniteDifferences and not our EnzymeRules
                # see comments in https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/589
                if !(Grid <: AbstractReducedGrid) & fd_tests[i_grid]
                    trunc = 5
                    spectrum = Spectrum(trunc, one_degree_more = true)
                    grid = Grid(SpeedyTransforms.get_nlat_half(trunc, grid_dealiasing[i_grid]))
                    S = SpectralTransform(spectrum, grid)
                    field = rand(grid)
                    f_north = S.scratch_memory.north
                    f_south = S.scratch_memory.south

                    # forward transform
                    test_reverse(
                        SpeedyTransforms._fourier!, Const,
                        (f_north, Duplicated), (f_south, Duplicated), (field, Duplicated), (S, Const);
                        fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                    )

                    # inverse transform
                    field = zero(field)
                    test_reverse(
                        SpeedyTransforms._fourier!, Const,
                        (field, Duplicated), (f_north, Duplicated), (f_south, Duplicated), (S, Const);
                        fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                    )
                end
            end
        end

        # `runtime_activity` is needed throughout the forward tests because EnzymeTestUtils' own
        # call wrapper stores constant memory into a differentiable tuple, which trips Enzyme's
        # static activity analysis (the reverse tests above hit the same thing, hence the
        # `set_runtime_activity(Reverse)` in the `autodiff` calls further down).
        @testset "EnzymeTestUtils forward rule test" begin
            for (i_grid, Grid) in enumerate(grid_types)

                # same reduced-grid caveat as the reverse tests above
                if !(Grid <: AbstractReducedGrid) & fd_tests[i_grid]
                    trunc = 5
                    spectrum = Spectrum(trunc, one_degree_more = true)
                    grid = Grid(SpeedyTransforms.get_nlat_half(trunc, grid_dealiasing[i_grid]))
                    S = SpectralTransform(spectrum, grid)
                    field = rand(grid)
                    f_north = zero(S.scratch_memory.north)
                    f_south = zero(S.scratch_memory.south)

                    # forward transform
                    test_forward(
                        fourier_analysis!, Duplicated,
                        (f_north, Duplicated), (f_south, Duplicated), (field, Duplicated), (S, Const);
                        fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                        runtime_activity = true,
                    )

                    # inverse transform
                    test_forward(
                        fourier_synthesis!, Duplicated,
                        (zero(field), Duplicated), (f_north, Duplicated), (f_south, Duplicated), (S, Const);
                        fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                        runtime_activity = true,
                    )
                end
            end
        end
    end
    @testset "chunked transform Enzyme rules" begin
        # SpeedyTransformsEnzymeExt defines custom reverse rules for the CHUNKED (unplanned-K)
        # CPU transform that apply the analytic adjoint of the (linear) spectral transform.
        # The chunked-path gradient must match the batched path, which native Enzyme + the
        # _fourier! rules already handle correctly. (Differentiating the chunk loop itself is
        # unsafe: Enzyme mis-builds the per-chunk view shadows — degenerate (0,1) shadow / GC
        # corruption on 1.10 — and reuses the last iteration's shadow for all chunks.)
        trunc = 5
        spectrum = Spectrum(trunc, one_degree_more = true)
        grid = FullGaussianGrid(SpeedyTransforms.get_nlat_half(trunc, grid_dealiasing[1]))
        NL = 4
        S_chunked = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1])
        S_batched = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1, NL])
        @test SpeedyTransforms._needs_chunking(NL, S_chunked)
        @test !SpeedyTransforms._needs_chunking(NL, S_batched)

        # EnzymeTestUtils reverse-rule check against finite differences, for the CHUNKED transform
        # (exercises the analytic-adjoint transform! rule; scratch Const so FD doesn't perturb it).
        # Wrapped to return `nothing` (like `_fourier!`) — otherwise ETU treats the returned mutated
        # output array as an extra differentiable output and its FD Jacobian goes out of bounds.
        transform_s2g!(field, coeffs, scratch, S) = (transform!(field, coeffs, scratch, S); nothing)
        transform_g2s!(coeffs, field, scratch, S) = (transform!(coeffs, field, scratch, S); nothing)
        @testset "EnzymeTestUtils reverse rule test" begin
            let field = zeros(Float32, grid, NL), coeffs = rand(ComplexF32, spectrum, NL)
                test_reverse(
                    transform_s2g!, Const,
                    (field, Duplicated), (coeffs, Duplicated),
                    (deepcopy(S_chunked.scratch_memory), Const), (S_chunked, Const);
                    fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                )
            end
            let coeffs = zeros(ComplexF32, spectrum, NL), field = rand(Float32, grid, NL)
                test_reverse(
                    transform_g2s!, Const,
                    (coeffs, Duplicated), (field, Duplicated),
                    (deepcopy(S_chunked.scratch_memory), Const), (S_chunked, Const);
                    fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                )
            end
        end

        # spec -> grid: vjp w.r.t coeffs must agree between chunked and batched transforms
        coeffs0 = rand(ComplexF32, spectrum, NL)
        dfield0 = rand(Float32, grid, NL)
        dcoeffs = map((S_chunked, S_batched)) do S
            coeffs = deepcopy(coeffs0)
            field = zeros(Float32, grid, NL)
            dfield = deepcopy(dfield0)
            dc = make_zero(coeffs)
            autodiff(
                set_runtime_activity(Reverse), transform!, Const,
                Duplicated(field, dfield), Duplicated(coeffs, dc),
                Duplicated(deepcopy(S.scratch_memory), make_zero(deepcopy(S.scratch_memory))), Const(S),
            )
            dc
        end
        @test all(isfinite, dcoeffs[1].data)
        @test any(!iszero, dcoeffs[1].data)
        @test isapprox(dcoeffs[1].data, dcoeffs[2].data, rtol = 1.0e-4)

        # the rule must ACCUMULATE into the (input) coeffs cotangent, not overwrite it: seeding a
        # pre-existing gradient `base_c` must yield base_c + pullback (guards the reset=false path)
        base_c = rand(ComplexF32, spectrum, NL)
        dc_acc = deepcopy(base_c)
        autodiff(
            set_runtime_activity(Reverse), transform!, Const,
            Duplicated(zeros(Float32, grid, NL), deepcopy(dfield0)), Duplicated(deepcopy(coeffs0), dc_acc),
            Duplicated(deepcopy(S_chunked.scratch_memory), make_zero(deepcopy(S_chunked.scratch_memory))), Const(S_chunked),
        )
        @test isapprox(dc_acc.data, base_c.data .+ dcoeffs[1].data, rtol = 1.0e-4)

        # grid -> spec: vjp w.r.t field must agree between chunked and batched transforms
        field0 = rand(Float32, grid, NL)
        dcoeffs0 = rand(ComplexF32, spectrum, NL)
        dfields = map((S_chunked, S_batched)) do S
            field = deepcopy(field0)
            coeffs = zeros(ComplexF32, spectrum, NL)
            dcoeffs_seed = deepcopy(dcoeffs0)
            df = make_zero(field)
            autodiff(
                set_runtime_activity(Reverse), transform!, Const,
                Duplicated(coeffs, dcoeffs_seed), Duplicated(field, df),
                Duplicated(deepcopy(S.scratch_memory), make_zero(deepcopy(S.scratch_memory))), Const(S),
            )
            df
        end
        @test all(isfinite, dfields[1].data)
        @test any(!iszero, dfields[1].data)
        @test isapprox(dfields[1].data, dfields[2].data, rtol = 1.0e-4)

        # accumulate into the (input) field cotangent, not overwrite (guards the add=true path)
        base_f = rand(Float32, grid, NL)
        df_acc = deepcopy(base_f)
        autodiff(
            set_runtime_activity(Reverse), transform!, Const,
            Duplicated(zeros(ComplexF32, spectrum, NL), deepcopy(dcoeffs0)), Duplicated(deepcopy(field0), df_acc),
            Duplicated(deepcopy(S_chunked.scratch_memory), make_zero(deepcopy(S_chunked.scratch_memory))), Const(S_chunked),
        )
        @test isapprox(df_acc.data, base_f.data .+ dfields[1].data, rtol = 1.0e-4)

        # FORWARD MODE. `transform!` already returns the array it mutates, as `test_forward`
        # requires, so no wrapper is needed here. Only the chunked transform is checked against
        # finite differences — the rule is one code path for chunked and batched alike — while
        # chunked/batched agreement is checked exactly in the next testset.
        @testset "EnzymeTestUtils forward rule test" begin
            let field = zeros(Float32, grid, NL), coeffs = rand(ComplexF32, spectrum, NL)
                test_forward(
                    transform!, Duplicated,
                    (field, Duplicated), (coeffs, Duplicated),
                    (deepcopy(S_chunked.scratch_memory), Const), (S_chunked, Const);
                    fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                    runtime_activity = true,
                )
            end
            let coeffs = zeros(ComplexF32, spectrum, NL), field = rand(Float32, grid, NL)
                test_forward(
                    transform!, Duplicated,
                    (coeffs, Duplicated), (field, Duplicated),
                    (deepcopy(S_chunked.scratch_memory), Const), (S_chunked, Const);
                    fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
                    runtime_activity = true,
                )
            end
        end

        # The transform is linear and `S` is inactive, so the jvp must be EXACTLY the transform of
        # the tangent. That is a much tighter check than finite differences and costs almost
        # nothing, so it is run for both the chunked and the batched path.
        @testset "forward rule jvp == transform(tangent)" begin
            for S in (S_chunked, S_batched)
                # spec -> grid
                coeffs, dcoeffs = rand(ComplexF32, spectrum, NL), rand(ComplexF32, spectrum, NL)
                field, dfield = zeros(Float32, grid, NL), zeros(Float32, grid, NL)
                autodiff(
                    Forward, transform!, Const,
                    Duplicated(field, dfield), Duplicated(coeffs, dcoeffs),
                    Const(deepcopy(S.scratch_memory)), Const(S),
                )
                @test any(!iszero, dfield)
                @test dfield ≈ transform!(zeros(Float32, grid, NL), dcoeffs, deepcopy(S.scratch_memory), S)
                @test field ≈ transform!(zeros(Float32, grid, NL), coeffs, deepcopy(S.scratch_memory), S)

                # grid -> spec
                field2, dfield2 = rand(Float32, grid, NL), rand(Float32, grid, NL)
                coeffs2, dcoeffs2 = zeros(ComplexF32, spectrum, NL), zeros(ComplexF32, spectrum, NL)
                autodiff(
                    Forward, transform!, Const,
                    Duplicated(coeffs2, dcoeffs2), Duplicated(field2, dfield2),
                    Const(deepcopy(S.scratch_memory)), Const(S),
                )
                @test any(!iszero, dcoeffs2)
                @test dcoeffs2 ≈ transform!(zeros(ComplexF32, spectrum, NL), dfield2, deepcopy(S.scratch_memory), S)
                @test coeffs2 ≈ transform!(zeros(ComplexF32, spectrum, NL), field2, deepcopy(S.scratch_memory), S)
            end

            # width > 1 (BatchDuplicated): every tangent must be transformed independently
            coeffs = rand(ComplexF32, spectrum, NL)
            dcoeffs = (rand(ComplexF32, spectrum, NL), rand(ComplexF32, spectrum, NL))
            refs = map(dc -> transform!(zeros(Float32, grid, NL), dc, deepcopy(S_batched.scratch_memory), S_batched), dcoeffs)
            field = zeros(Float32, grid, NL)
            dfields = (zeros(Float32, grid, NL), zeros(Float32, grid, NL))
            autodiff(
                Forward, transform!, Const,
                BatchDuplicated(field, dfields), BatchDuplicated(coeffs, dcoeffs),
                Const(deepcopy(S_batched.scratch_memory)), Const(S_batched),
            )
            @test dfields[1] ≈ refs[1]
            @test dfields[2] ≈ refs[2]

            # an inactive (Const) input carries a zero tangent, so the output tangent is zeroed
            dfield = zeros(Float32, grid, NL)
            dfield .= 7     # pre-seeded, must not survive
            autodiff(
                Forward, transform!, Const,
                Duplicated(zeros(Float32, grid, NL), dfield), Const(rand(ComplexF32, spectrum, NL)),
                Const(deepcopy(S_batched.scratch_memory)), Const(S_batched),
            )
            @test all(iszero, dfield)
        end
    end
    # The reverse rules end with `make_zero!` on their output cotangent, which is correct for a
    # plain array: the primal overwrites that argument, so its input-side cotangent is zero (the
    # `test_reverse` checks above pin exactly that against finite differences).
    #
    # In SpeedyWeather the transform arguments are NOT plain arrays — they are views into a shared
    # "fused" parent buffer that several variables live in. Zeroing then reaches beyond the slice
    # this call overwrote. These tests exercise that layout directly, which the plain-array tests
    # above cannot: they check that a cotangent parked in a DISJOINT slice of the same parent
    # survives a transform on a neighbouring slice, and that the transform's own slice still
    # pulls back correctly.
    @testset "reverse rules on view-backed (fused) arrays" begin
        trunc = 5
        spectrum = Spectrum(trunc, one_degree_more = true)
        grid = FullGaussianGrid(SpeedyTransforms.get_nlat_half(trunc, grid_dealiasing[1]))
        NL = 2
        S = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1, NL])

        # one parent buffer holding two independent layer-slices, as the fused variables do
        spec_parent = rand(ComplexF32, spectrum, 2NL)
        grid_parent = rand(Float32, grid, 2NL)
        spec_a = wrapped_view(spec_parent, :, 1:NL)         # the transform's argument
        spec_b = wrapped_view(spec_parent, :, (NL + 1):2NL) # a neighbour that must be left alone
        field_a = wrapped_view(grid_parent, :, 1:NL)
        field_b = wrapped_view(grid_parent, :, (NL + 1):2NL)

        # spec -> grid: cotangent on the neighbouring grid slice must survive
        let dspec = make_zero(spec_parent), dgrid = make_zero(grid_parent)
            dfield_b = wrapped_view(dgrid, :, (NL + 1):2NL)
            fill!(dfield_b.data, 1)                        # mark the neighbour's cotangent
            marker = copy(dfield_b.data)
            autodiff(
                set_runtime_activity(Reverse), transform!, Const,
                Duplicated(field_a, wrapped_view(dgrid, :, 1:NL)),
                Duplicated(spec_a, wrapped_view(dspec, :, 1:NL)),
                Const(deepcopy(S.scratch_memory)), Const(S),
            )
            @test wrapped_view(dgrid, :, (NL + 1):2NL).data == marker   # neighbour untouched
            @test any(!iszero, wrapped_view(dspec, :, 1:NL).data)       # own slice pulled back
        end

        # grid -> spec: same, with the roles swapped
        let dspec = make_zero(spec_parent), dgrid = make_zero(grid_parent)
            dspec_b = wrapped_view(dspec, :, (NL + 1):2NL)
            fill!(dspec_b.data, 1 + im)
            marker = copy(dspec_b.data)
            autodiff(
                set_runtime_activity(Reverse), transform!, Const,
                Duplicated(spec_a, wrapped_view(dspec, :, 1:NL)),
                Duplicated(field_a, wrapped_view(dgrid, :, 1:NL)),
                Const(deepcopy(S.scratch_memory)), Const(S),
            )
            @test wrapped_view(dspec, :, (NL + 1):2NL).data == marker   # neighbour untouched
            @test any(!iszero, wrapped_view(dgrid, :, 1:NL).data)       # own slice pulled back
        end
    end

    @testset "Complete Transform ChainRules" begin
        # WIP
    end
end
