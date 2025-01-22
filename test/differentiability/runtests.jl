# Here we test the differentiability of the model by comparing Enzyme's results to finite difference differenation 
# These tests are relatively expensive, and also not strictly necessary to perform at every commit, so they 
# are only part of the extended test set 
using SpeedyWeather, EnzymeTestUtils, Enzyme, FiniteDifferences, Test
import EnzymeTestUtils: test_approx
import FiniteDifferences: j′vp, grad, central_fdm
import AbstractFFTs

grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms 
grid_dealiasing = [2, 3]
fd_tests = [true, true] 

# UTILITIES 
include("test_utils.jl")
include("timestep_utils.jl")

# TESTS 
include("speedy_transforms.jl")
include("barotropic.jl")