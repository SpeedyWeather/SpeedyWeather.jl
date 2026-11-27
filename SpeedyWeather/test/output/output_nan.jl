using NCDatasets

# NaN bit patterns that the mantissa rounding of BitInformation.round! turns into ±0 or ±Inf:
# all trailing payload bits set (0x7fffffff is the NaN that CUDA arithmetic produces) or only
# low payload bits set, besides the default NaNs of Julia (0x7fc00000) and x86 (0xffc00000)
const NAN_BIT_PATTERNS = (
    Float16 => (0x7e00, 0xfe00, 0x7fff, 0xffff, 0x7c01, 0x7c1f),
    Float32 => (0x7fc00000, 0xffc00000, 0x7fffffff, 0xffffffff, 0x7f800001, 0x7f8003ff),
    Float64 => (0x7ff8000000000000, 0xfff8000000000000, 0x7fffffffffffffff, 0xffffffffffffffff, 0x7ff00000000003ff),
)

# output variable that writes `field` as it is, to inject NaN bit patterns into the output
@kwdef mutable struct NaNBitPatternOutput{F} <: SpeedyWeather.AbstractOutputVariable
    name::String = "nan_patterns"
    unit::String = "1"
    long_name::String = "ones with NaNs of different bit patterns"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    keepbits::Int = 7
    field::F
end

SpeedyWeather.OutputVariables.path(variable::NaNBitPatternOutput, simulation) = variable.field

@testset "Mantissa rounding keeps NaNs" begin
    for (T, patterns) in NAN_BIT_PATTERNS
        U = Base.uinttype(T)
        values = T[1.2345, -6.789e-3, 0, -0.0, 3, Inf, -Inf, floatmax(T)]
        nans = [reinterpret(T, U(pattern)) for pattern in patterns]
        @test all(isnan, nans)

        for keepbits in (0, 1, 7, Base.significand_bits(T) - 1, Base.significand_bits(T))
            x = vcat(values, nans)
            SpeedyWeather.round_mantissa!(x, keepbits)

            # every NaN stays NaN, written as Julia's NaN
            @test all(xᵢ -> reinterpret(U, xᵢ) == reinterpret(U, T(NaN)), x[(length(values) + 1):end])

            # all other values as BitInformation's rounding
            @test all(reinterpret.(U, x[1:length(values)]) .== reinterpret.(U, round.(values, keepbits)))
        end
    end
end

@testset "NetCDF output keeps NaNs of any bit pattern" begin
    path = mktempdir(pwd(), prefix = "tmp_output_nan_")
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)

    # ones with every 10th point a NaN of CUDA's bit pattern
    field = ones(Float32, spectral_grid.grid)
    cuda_nan = reinterpret(Float32, 0x7fffffff)
    field.data[1:10:end] .= cuda_nan

    output = NetCDFOutput(spectral_grid; path)
    add!(output, NaNBitPatternOutput(; field))
    model = BarotropicModel(spectral_grid; output)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, steps = 1, output = true)

    written = NCDataset(joinpath(output.run_path, output.filename)) do ds
        ds["nan_patterns"].var[:, :, end]
    end

    # the interpolation onto the output grid carries the NaNs (and their bit pattern)
    # to the neighbouring points, every other point is 1 (no NaN written as ±0 or ±Inf)
    @test 0 < count(isnan, written) < length(written)
    @test all(x -> isnan(x) || x ≈ 1, written)
end
