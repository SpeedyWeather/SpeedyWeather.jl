# Unit tests for the custom EnzymeRules defined in SpeedyWeatherEnzymeExt.
#
# `get_step(var, step)` returns a view selecting a leapfrog/time step along the last dimension.
# It has a custom reverse rule (see ext/SpeedyWeatherEnzymeExt.jl) so the runtime-offset view
# does not trip Enzyme's type analysis on Julia ≥ 1.11.
#
# `get_step` returns a view that *aliases* its input, which `test_reverse` cannot seed directly
# (the return shadow would be the same memory as the checked input shadow). We therefore exercise
# the rule through two non-aliasing wrappers — `copy` (materialise the step view) and a scalar
# reduction — and `test_reverse` checks the rule's vjp against finite differences for every step
# index and array kind.
@testset "Differentiability: get_step custom EnzymeRule" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 3)
    spectrum = spectral_grid.spectrum
    grid = spectral_grid.grid

    step_copy(var, step) = copy(SpeedyWeather.get_step(var, step))   # materialised step → non-aliasing array return
    step_sum(var, step) = sum(abs2, SpeedyWeather.get_step(var, step))  # scalar (Active) return

    @testset "3D LowerTriangularArray (coeffs × nlayers × steps)" begin
        var = rand(ComplexF32, spectrum, 3, 2)
        for step in 1:2
            test_reverse(step_copy, Duplicated, (var, Duplicated), (step, Const))
            test_reverse(step_sum, Active, (var, Duplicated), (step, Const))
        end
    end

    @testset "2D LowerTriangularArray (coeffs × steps)" begin
        var = rand(ComplexF32, spectrum, 2)
        for step in 1:2
            test_reverse(step_copy, Duplicated, (var, Duplicated), (step, Const))
            test_reverse(step_sum, Active, (var, Duplicated), (step, Const))
        end
    end

    @testset "3D Field (grid × nlayers × steps)" begin
        var = rand(Float32, grid, 3, 2)
        for step in 1:2
            test_reverse(step_copy, Duplicated, (var, Duplicated), (step, Const))
            test_reverse(step_sum, Active, (var, Duplicated), (step, Const))
        end
    end

    @testset "2D Field (grid × steps)" begin
        var = rand(Float32, grid, 2)
        for step in 1:2
            test_reverse(step_copy, Duplicated, (var, Duplicated), (step, Const))
            test_reverse(step_sum, Active, (var, Duplicated), (step, Const))
        end
    end
end
