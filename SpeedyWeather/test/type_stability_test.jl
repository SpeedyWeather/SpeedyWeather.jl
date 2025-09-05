# Test type stability of SpectralTransform
# Run this from the package directory with: julia --project=. test/type_stability_test.jl

import Pkg 
Pkg.activate(".")
using Test
using InteractiveUtils  # For @code_warntype
using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using SpeedyWeather.LowerTriangularArrays

# Create a simple test function that uses SpectralTransform
function test_spectral_transform_stability()
    # Create a grid and spectrum
    grid = FullGaussianGrid(32)  # 32 latitude rings per hemisphere
    spectrum = Spectrum(42, 42)  # T41 truncation (42x42 coefficients, 0-based indexing)
    
    # Create a SpectralTransform
    S = SpectralTransform(spectrum, grid)
    
    # Create a field and coefficients
    field = zeros(Float64, grid)
    coeffs = zeros(Complex{Float64}, spectrum)
    
    # Set some values in the field
    lons, lats = RingGrids.get_lonlats(field.grid)
    for i in eachindex(field)
        field[i] = sin(lons[i]) * cos(lats[i])
    end
    
    # Perform a transform
    transform!(coeffs, field, S)
    
    # Transform back
    field_reconstructed = transform(coeffs, S)
    
    return field_reconstructed
end

# Test the function with @code_warntype
println("Testing SpectralTransform type stability...")
println("\nTesting transform! (grid to spectral):")
@code_warntype transform!(zeros(Complex{Float64}, Spectrum(42, 42)), 
                         zeros(Float64, FullGaussianGrid(32)), 
                         SpectralTransform(Spectrum(42, 42), FullGaussianGrid(32)))

println("\nTesting transform (spectral to grid):")
@code_warntype transform(zeros(Complex{Float64}, Spectrum(42, 42)), 
                        SpectralTransform(Spectrum(42, 42), FullGaussianGrid(32)))

# Also test the core transform functions
println("\nTesting _fourier! (grid to spectral):")
function test_fourier_stability()
    grid = FullGaussianGrid(32)
    spectrum = Spectrum(42, 42)
    S = SpectralTransform(spectrum, grid)
    field = zeros(Float64, grid)
    f_north = S.scratch_memory_north
    f_south = S.scratch_memory_south
    SpeedyTransforms._fourier!(f_north, f_south, field, S)
end
@code_warntype test_fourier_stability()

println("\nTesting _legendre! (grid to spectral):")
function test_legendre_stability()
    grid = FullGaussianGrid(32)
    spectrum = Spectrum(42, 42)
    S = SpectralTransform(spectrum, grid)
    coeffs = zeros(Complex{Float64}, spectrum)
    f_north = S.scratch_memory_north
    f_south = S.scratch_memory_south
    SpeedyTransforms._legendre!(coeffs, f_north, f_south, S)
end
@code_warntype test_legendre_stability()

# Run the test function to see if it works correctly
test_spectral_transform_stability()
println("Test completed successfully!")
