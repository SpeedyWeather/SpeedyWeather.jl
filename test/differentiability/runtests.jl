# Here we test the differentiability of the model by comparing Enzyme's results to finite differences. 
# These tests are relatively expensive, and also not strictly necessary to perform at every commit, so they 
# are only part of the extended test set 
using SpeedyWeather
using ComponentArrays
using EnzymeTestUtils, Enzyme, FiniteDifferences, StatsBase, Test

import EnzymeTestUtils: test_approx
import FiniteDifferences: j′vp, grad, central_fdm
import AbstractFFTs

# UTILITIES 
include("test_utils.jl")
include("timestep_utils.jl")

# TESTS
include("speedy_transforms.jl")
include("barotropic.jl")
#include("primitivewet.jl")
include("timestepping.jl")