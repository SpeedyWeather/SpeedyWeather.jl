# Test type stability of PrognosticVariables initialization
# Run this from the package directory with: julia --project=. test/prognostic_variables_stability_test.jl

import Pkg
Pkg.activate(".")
using Test
using InteractiveUtils  # For @code_warntype
using SpeedyWeather
using SpeedyWeather.SpeedyInternals.Architectures
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using SpeedyWeather.LowerTriangularArrays

# Create a simple test function that initializes PrognosticVariables
function test_prognostic_variables_stability()
    # Create a spectral grid with standard parameters
    spectral_grid = SpectralGrid(
        trunc=31,                    # T31 spectral truncation
        nlayers=8,                   # 8 vertical levels
    )

    # Initialize PrognosticVariables
    prognostic = PrognosticVariables(spectral_grid)
    
    return prognostic
end


# Test the function with @code_warntype
println("Testing PrognosticVariables initialization type stability...")
println("\nTesting SpectralGrid constructor:")

@code_warntype SpectralGrid(trunc=20)
@code_warntype test_prognostic_variables_stability()

# Run the test functions to see if they work correctly
println("\nRunning the test functions...")
test_prognostic_variables_stability()
println("Test completed successfully!")
