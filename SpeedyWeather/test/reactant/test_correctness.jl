"""
Test script comparing SpeedyWeather models on CPU vs Reactant.
Verifies that the Reactant implementation yields the same results as the CPU version.
"""

import Pkg
Pkg.activate(@__DIR__)

using SpeedyWeather
using Reactant
using Test
using LinearAlgebra
using CUDA
using Statistics: mean
import SpeedyWeather: ReactantDevice, first_timesteps!, later_timestep!

# ============================================================================
# Configuration
# ============================================================================
const TRUNC = 31            # spectral truncation
const NSTEPS = 10           # number of time steps to compare
const RTOL = 1.0e-3           # relative tolerance for comparison
const ATOL = 1.0e-8          # absolute tolerance for comparison

# ============================================================================
# Helper functions
# ============================================================================

"""Get nlayers for a given model type."""
nlayers_for_model(::Type{BarotropicModel}) = 1
nlayers_for_model(::Type{ShallowWaterModel}) = 1
nlayers_for_model(::Type{PrimitiveDryModel}) = 8
nlayers_for_model(::Type{PrimitiveWetModel}) = 8

"""Get default initial conditions for a given model type."""
function default_initial_conditions(::Type{BarotropicModel})
    return InitialConditions(; vordiv = ZeroInitially())
end

function default_initial_conditions(::Type{ShallowWaterModel})
    return InitialConditions(; vordiv = ZonalJet())
end

function default_initial_conditions(::Type{PrimitiveDryModel})
    return InitialConditions(; vordiv = ZonalWind())
end

function default_initial_conditions(::Type{PrimitiveWetModel})
    return InitialConditions(; vordiv = ZonalWind())
end

"""Create a CPU model of the given type."""
function create_cpu_model(ModelType::Type; trunc = TRUNC)
    nlayers = nlayers_for_model(ModelType)
    spectral_grid = SpectralGrid(; nlayers, trunc)
    initial_conditions = default_initial_conditions(ModelType)
    return ModelType(spectral_grid; initial_conditions)
end

"""Create a Reactant model of the given type."""
function create_reactant_model(ModelType::Type; trunc = TRUNC)
    nlayers = nlayers_for_model(ModelType)
    arch = SpeedyWeather.ReactantDevice()
    spectral_grid = SpectralGrid(; architecture = arch, nlayers, trunc)
    M = MatrixSpectralTransform(spectral_grid)
    initial_conditions = default_initial_conditions(ModelType)
    return ModelType(spectral_grid; spectral_transform = M, feedback = nothing, initial_conditions)
end

"""Synchronize variables from Reactant to CPU simulation."""
function sync_variables!(sim_cpu, sim_reactant)
    progn_cpu, _, _ = SpeedyWeather.unpack(sim_cpu)
    progn_reactant, _, _ = SpeedyWeather.unpack(sim_reactant)

    # Sync vorticity from Reactant to CPU (copy underlying data)
    copyto!(progn_cpu.vor.data, Array(progn_reactant.vor.data))

    #TODO: for the other models add more varibles
    return nothing
end

"""Compare prognostic variables between CPU and Reactant simulations."""
function compare_prognostic_variables(sim_cpu, sim_reactant; rtol = RTOL, atol = ATOL)
    progn_cpu, _, _ = SpeedyWeather.unpack(sim_cpu)
    progn_reactant, _, _ = SpeedyWeather.unpack(sim_reactant)

    results = Dict{Symbol, NamedTuple}()

    # Compare vorticity
    vor_cpu = Array(progn_cpu.vor[:, :, 2])
    vor_reactant = Array(progn_reactant.vor[:, :, 2])
    abs_diff = abs.(vor_cpu .- vor_reactant)
    results[:vor] = (
        max_abs_diff = maximum(abs_diff),
        mean_abs_diff = mean(abs_diff),
        matches = isapprox(vor_cpu, vor_reactant, rtol = rtol, atol = atol),
    )

    return results
end

"""Compare grid variables between CPU and Reactant simulations."""
function compare_grid_variables(sim_cpu, sim_reactant; rtol = RTOL, atol = ATOL)
    _, diagn_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    _, diagn_reactant, _ = SpeedyWeather.unpack(sim_reactant)

    results = Dict{Symbol, Bool}()

    # Compare u_grid
    u_cpu = Array(diagn_cpu.grid.u_grid)
    u_reactant = Array(diagn_reactant.grid.u_grid)
    results[:u_grid] = isapprox(u_cpu, u_reactant, rtol = rtol, atol = atol)

    # Compare v_grid
    v_cpu = Array(diagn_cpu.grid.v_grid)
    v_reactant = Array(diagn_reactant.grid.v_grid)
    results[:v_grid] = isapprox(v_cpu, v_reactant, rtol = rtol, atol = atol)

    # Compare vor_grid
    vor_cpu = Array(diagn_cpu.grid.vor_grid)
    vor_reactant = Array(diagn_reactant.grid.vor_grid)
    results[:vor_grid] = isapprox(vor_cpu, vor_reactant, rtol = rtol, atol = atol)

    return results
end

"""Compare tendencies between CPU and Reactant simulations after a single timestep."""
function compare_tendencies(sim_cpu, sim_reactant; rtol = RTOL, atol = ATOL)
    _, diagn_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    _, diagn_reactant, _ = SpeedyWeather.unpack(sim_reactant)

    results = Dict{Symbol, NamedTuple}()

    # Compare vorticity tendency
    vor_tend_cpu = Array(diagn_cpu.tendencies.vor_tend)
    vor_tend_reactant = Array(diagn_reactant.tendencies.vor_tend)
    abs_diff = abs.(vor_tend_cpu .- vor_tend_reactant)
    results[:vor_tend] = (
        max_abs_diff = maximum(abs_diff),
        mean_abs_diff = mean(abs_diff),
        matches = isapprox(vor_tend_cpu, vor_tend_reactant, rtol = rtol, atol = atol),
    )

    return results
end

"""Test tendencies after a single timestep on already-initialized simulations."""
function test_tendencies!(sim_cpu, sim_reactant, model_name; rtol = RTOL, atol = ATOL)
    println("\n" * "-"^60)
    println("Testing tendencies (single timestep)")
    println("-"^60)

    initialize!(sim_cpu, steps=1)
    initialize!(sim_reactant, steps=1)

    sync_variables!(sim_cpu, sim_reactant)

    # Run a single timestep to compute tendencies
    println("  Running CPU model...")
    SpeedyWeather.timestep!(sim_cpu)
    println("  Running Reactant model...")
    @jit SpeedyWeather.timestep!(sim_reactant)
    println("  ✓ Tendencies computed")

    # Compare tendencies
    tend_results = compare_tendencies(sim_cpu, sim_reactant; rtol, atol)

    println("\nVorticity tendency comparison:")
    println("  Max absolute difference:  $(tend_results[:vor_tend].max_abs_diff)")
    println("  Mean absolute difference: $(tend_results[:vor_tend].mean_abs_diff)")

    @testset "$model_name Tendency Comparison" begin
        @test tend_results[:vor_tend].matches
    end

    return tend_results
end

"""Test prognostic and grid variables after running for nsteps on already-initialized simulations."""
function test_time_stepping!(sim_cpu, sim_reactant, model_name; nsteps = NSTEPS, rtol = RTOL, atol = ATOL)
    println("\n" * "-"^60)
    println("Testing time stepping ($nsteps steps)")
    println("-"^60)
    
    sync_variables!(sim_cpu, sim_reactant)

    # Run time stepping
    println("  Running CPU model...")
    SpeedyWeather.run!(sim_cpu; steps = nsteps)

    println("  Running Reactant model...")
    SpeedyWeather.run!(sim_reactant; steps = nsteps)
    println("  ✓ Time stepping completed")

    # Compare results
    progn_results = compare_prognostic_variables(sim_cpu, sim_reactant; rtol, atol)
    grid_results = compare_grid_variables(sim_cpu, sim_reactant; rtol, atol)

    println("\nVorticity comparison after $nsteps steps:")
    println("  Max absolute difference:  $(progn_results[:vor].max_abs_diff)")
    println("  Mean absolute difference: $(progn_results[:vor].mean_abs_diff)")

    @testset "$model_name Time Stepping" begin
        @testset "Prognostic variables" begin
            @test progn_results[:vor].matches
        end

        @testset "Grid variables" begin
            @test grid_results[:u_grid]
            @test grid_results[:v_grid]
            @test grid_results[:vor_grid]
        end
    end

    return (progn_results, grid_results)
end

"""Run all correctness tests for a given model type."""
function test_model(ModelType::Type; trunc = TRUNC, nsteps = NSTEPS, rtol = RTOL, atol = ATOL)
    model_name = string(ModelType)

    println("="^60)
    println("$model_name: CPU vs Reactant Correctness Tests")
    println("="^60)

    # Setup CPU model
    println("\n[1/3] Setting up CPU model...")
    model_cpu = create_cpu_model(ModelType; trunc)
    simulation_cpu = initialize!(model_cpu)
    println("  ✓ CPU model initialized (T$trunc)")

    # Setup Reactant model
    println("\n[2/3] Setting up Reactant model...")
    model_reactant = create_reactant_model(ModelType; trunc)
    simulation_reactant = initialize!(model_reactant)
    println("  ✓ Reactant model initialized")

    # spin up models a bit 
    println("\n[3/3] Spinning up models...")
    run!(simulation_cpu; period = Day(20))
    run!(simulation_reactant; period = Day(20))
    println("  ✓ Models spun up")

    # Run tests
    @testset "$model_name CPU vs Reactant" begin
        tend_results = test_tendencies!(simulation_cpu, simulation_reactant,     model_name; rtol, atol)
        stepping_results = test_time_stepping!(simulation_cpu, simulation_reactant, model_name; nsteps, rtol, atol)
    end

    println("\n" * "="^60)
    println("$model_name tests completed!")
    println("="^60)

    return nothing
end

# ============================================================================
# Run tests
# ============================================================================

test_model(BarotropicModel)
