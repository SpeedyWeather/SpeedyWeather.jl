# Driver: replicate runtests.jl preamble, then run only primitivewet (dyn-core group).
# Validates the looseTypeAnalysis Enzyme fix on the richer PrimitiveWet dynamics (Julia > 1.10).
using SpeedyWeather
using ComponentArrays
using EnzymeTestUtils, Enzyme, FiniteDifferences, StatsBase, Test

import EnzymeTestUtils: test_approx
import FiniteDifferences: j′vp, grad, central_fdm
import AbstractFFTs

const HERE = @__DIR__

include(joinpath(HERE, "test_utils.jl"))
include(joinpath(HERE, "timestep_utils.jl"))

println("=== Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme), " ===")

@testset "primitivewet" begin
    include(joinpath(HERE, "primitivewet.jl"))
end
