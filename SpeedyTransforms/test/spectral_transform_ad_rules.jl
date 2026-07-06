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
    end
    @testset "Complete Transform ChainRules" begin
        # WIP
    end
end
