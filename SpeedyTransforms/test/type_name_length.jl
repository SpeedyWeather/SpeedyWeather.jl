# Regression guard on the length of the `SpectralTransform` type name.
#
# Under NESTED Enzyme autodiff (reverse-over-forward / forward-over-forward) Enzyme names an
# argument value `arg.Const{SpectralTransform{…}}` after the FULL Julia type. LLVM caps local value
# names at ~1024 characters and truncates beyond that; two identically-typed `Const{SpectralTransform}`
# values (one per nesting level) then truncate to the same name and collide with
#
#     LLVM error: multiple definition of local value named 'arg.Const{SpectralTransform{…}}'
#
# which blocks nested AD entirely. Testing the actual nested AD in the unit tests takes a long time, so we test here 
# just the length of the type name.
@testset "SpectralTransform type name length" begin
    ENZYME_VALUE_NAME_CAP = 1024      # LLVM local value-name truncation limit
    PREFIX = length("arg.Const{") + length("}")   # what Enzyme wraps the type name in

    for Grid in (FullGaussianGrid, OctahedralGaussianGrid), NF in (Float32, Float64)
        grid = Grid(SpeedyTransforms.get_nlat_half(31, 2))
        S = SpectralTransform(Spectrum(31), grid; NF)
        value_name_length = length(string(typeof(S))) + PREFIX
        @test value_name_length < ENZYME_VALUE_NAME_CAP
        # surface how close we are, so a shrinking margin is visible before it becomes a failure
        @info "SpectralTransform value-name length" Grid NF value_name_length margin = ENZYME_VALUE_NAME_CAP - value_name_length
    end
end
