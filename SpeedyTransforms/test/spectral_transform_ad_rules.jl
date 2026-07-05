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
    @testset "wrapped_view Enzyme rules" begin
        # alias rules in SpeedyTransformsEnzymeExt: the shadow of a wrapped_view is built
        # explicitly as the same view of the parent's shadow (fixes the reverse pass of
        # chunked/fused transforms where Enzyme mis-constructs that shadow)
        trunc = 5
        spectrum = Spectrum(trunc, one_degree_more = true)
        grid = FullGaussianGrid(SpeedyTransforms.get_nlat_half(trunc, grid_dealiasing[1]))

        # mutating consumers routing data through wrapped_view on both wrapper types;
        # gradients must match finite differences
        function field_chunk_double!(out, in)
            a = SpeedyTransforms.wrapped_view(in, :, 3:4)
            b = SpeedyTransforms.wrapped_view(out, :, 1:2)
            b.data .= 2 .* a.data
            return nothing
        end
        field_in = rand(Float32, grid, 4)
        field_out = zeros(Float32, grid, 4)
        test_reverse(
            field_chunk_double!, Const, (field_out, Duplicated), (field_in, Duplicated);
            fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
        )

        function lta_chunk_double!(out, in)
            a = SpeedyTransforms.wrapped_view(in, :, 3:4)
            b = SpeedyTransforms.wrapped_view(out, :, 1:2)
            b.data .= 2 .* a.data
            return nothing
        end
        lta_in = rand(Float32, spectrum, 4)
        lta_out = zeros(Float32, spectrum, 4)
        test_reverse(
            lta_chunk_double!, Const, (lta_out, Duplicated), (lta_in, Duplicated);
            fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
        )

        # integration: the CHUNKED transform (unplanned K routes through wrapped_view chunks)
        # must produce the same gradient as the batched transform (no chunking)
        NL = 4
        S_chunked = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1])
        S_batched = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1, NL])
        @test SpeedyTransforms._needs_chunking(NL, S_chunked)

        dcoeffs = map((S_chunked, S_batched)) do S
            coeffs = rand(ComplexF32, spectrum, NL)
            field = zeros(Float32, grid, NL)
            scratch = deepcopy(S.scratch_memory)
            dfield = zero(field)
            dfield .= 1
            dc = make_zero(coeffs)
            autodiff(
                set_runtime_activity(Reverse), transform!, Const,
                Duplicated(field, dfield), Duplicated(coeffs, dc),
                Duplicated(scratch, make_zero(scratch)), Const(S),
            )
            dc
        end
        @test all(isfinite, dcoeffs[1].data)
        @test any(!iszero, dcoeffs[1].data)
        @test isapprox(dcoeffs[1].data, dcoeffs[2].data, rtol = 1.0e-4)
    end
    @testset "Complete Transform ChainRules" begin
        # WIP
    end
end
