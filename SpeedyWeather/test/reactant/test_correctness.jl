"""
Test script comparing Barotropic model on CPU vs Reactant.
Verifies that the Reactant implementation yields the same results as the CPU version.
"""

import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path = joinpath(@__DIR__, "..", ".."))

using SpeedyWeather
using Reactant
using Test
using LinearAlgebra
using Statistics: mean

# Configuration
const TRUNC = 31            # spectral truncation
const NLAYERS = 1           # barotropic = 1 layer
const NSTEPS = 10           # number of time steps to compare
const RTOL = 1e-5           # relative tolerance for comparison
const ATOL = 1e-10          # absolute tolerance for comparison

println("=" ^ 60)
println("Barotropic Model: CPU vs Reactant Comparison Test")
println("=" ^ 60)

# ============================================================================
# Setup CPU Model
# ============================================================================
println("\n[1/4] Setting up CPU model...")

spectral_grid_cpu = SpectralGrid(; nlayers = NLAYERS, trunc = TRUNC)

# Use ZonalJet for deterministic initial conditions
model_cpu = BarotropicModel(
    spectral_grid_cpu;
    initial_conditions = InitialConditions(; vordiv = ZonalJet()),
)

simulation_cpu = initialize!(model_cpu)
initialize!(simulation_cpu; steps = NSTEPS, output = false)

println("  ✓ CPU model initialized")
println("  Resolution: T$(TRUNC), $(spectral_grid_cpu.nlat) latitudes")

# ============================================================================
# Setup Reactant Model
# ============================================================================
println("\n[2/4] Setting up Reactant model...")

arch_reactant = SpeedyWeather.ReactantDevice()
spectral_grid_reactant = SpectralGrid(; architecture = arch_reactant, nlayers = NLAYERS, trunc = TRUNC)

# Use MatrixSpectralTransform for Reactant compatibility
M = MatrixSpectralTransform(spectral_grid_reactant)

model_reactant = BarotropicModel(
    spectral_grid_reactant;
    spectral_transform = M,
    feedback = nothing,  # disable feedback for Reactant
    initial_conditions = InitialConditions(; vordiv = ZonalJet()),
)

simulation_reactant = initialize!(model_reactant)
initialize!(simulation_reactant; steps = NSTEPS, output = false)

println("  ✓ Reactant model initialized")

# ============================================================================
# Copy CPU initial state to Reactant model for fair comparison
# ============================================================================
println("\n[3/4] Synchronizing initial conditions...")

progn_cpu, diagn_cpu, _ = SpeedyWeather.unpack(simulation_cpu)
progn_reactant, diagn_reactant, _ = SpeedyWeather.unpack(simulation_reactant)

# Copy vorticity from CPU to Reactant (both leapfrog steps)
vor_cpu_1 = Array(progn_cpu.vor[:, :, 1])
vor_cpu_2 = Array(progn_cpu.vor[:, :, 2])

# Use Reactant.to_rarray to transfer data
progn_reactant.vor[:, :, 1] .= Reactant.to_rarray(vor_cpu_1)
progn_reactant.vor[:, :, 2] .= Reactant.to_rarray(vor_cpu_2)

println("  ✓ Initial vorticity synchronized")

# Verify initial state matches
vor_reactant_check = Array(progn_reactant.vor[:, :, 1])
@assert isapprox(vor_cpu_1, vor_reactant_check, rtol = 1e-10) "Initial state mismatch!"
println("  ✓ Initial state verified")

# ============================================================================
# Run time stepping and compare
# ============================================================================
println("\n[4/4] Running time stepping comparison...")

# Store initial state for reference
vor_initial = copy(vor_cpu_1)

# Run CPU model
println("  Running CPU model for $NSTEPS steps...")
SpeedyWeather.time_stepping!(simulation_cpu)

# Run Reactant model
println("  Running Reactant model for $NSTEPS steps...")
SpeedyWeather.time_stepping!(simulation_reactant)

# ============================================================================
# Compare Results
# ============================================================================
println("\n" * "=" ^ 60)
println("RESULTS")
println("=" ^ 60)

# Extract final vorticity
vor_cpu_final = Array(progn_cpu.vor[:, :, 2])
vor_reactant_final = Array(progn_reactant.vor[:, :, 2])

# Compute differences
abs_diff = abs.(vor_cpu_final .- vor_reactant_final)
max_abs_diff = maximum(abs_diff)
mean_abs_diff = mean(abs_diff)

# Relative difference (avoid division by zero)
rel_diff = abs_diff ./ (abs.(vor_cpu_final) .+ ATOL)
max_rel_diff = maximum(rel_diff)
mean_rel_diff = mean(rel_diff)

println("\nVorticity comparison after $NSTEPS steps:")
println("  Max absolute difference:  $(max_abs_diff)")
println("  Mean absolute difference: $(mean_abs_diff)")
println("  Max relative difference:  $(max_rel_diff)")
println("  Mean relative difference: $(mean_rel_diff)")

# Test assertions
@testset "Barotropic CPU vs Reactant" begin
    @testset "Vorticity" begin
        @test isapprox(vor_cpu_final, vor_reactant_final, rtol = RTOL, atol = ATOL)
        @test max_abs_diff < RTOL
    end

    # Also compare grid variables if available
    @testset "Grid variables" begin
        u_cpu = Array(diagn_cpu.grid.u_grid)
        u_reactant = Array(diagn_reactant.grid.u_grid)
        @test isapprox(u_cpu, u_reactant, rtol = RTOL, atol = ATOL)

        v_cpu = Array(diagn_cpu.grid.v_grid)
        v_reactant = Array(diagn_reactant.grid.v_grid)
        @test isapprox(v_cpu, v_reactant, rtol = RTOL, atol = ATOL)

        vor_grid_cpu = Array(diagn_cpu.grid.vor_grid)
        vor_grid_reactant = Array(diagn_reactant.grid.vor_grid)
        @test isapprox(vor_grid_cpu, vor_grid_reactant, rtol = RTOL, atol = ATOL)
    end
end

println("\n" * "=" ^ 60)
println("Test completed!")
println("=" ^ 60)
