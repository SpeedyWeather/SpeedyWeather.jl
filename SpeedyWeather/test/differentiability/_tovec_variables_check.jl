# Validate the custom to_vec(::Variables): round-trip, perturbation propagates through fused views,
# and a small FiniteDifferences j′vp runs without the Array->SubArray error.
# Run: julia --project=SpeedyWeather/test/differentiability <file>
using SpeedyWeather
using FiniteDifferences, Test
import FiniteDifferences: to_vec, j′vp, central_fdm

spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
sim = initialize!(model)
vars = sim.variables

# count view leaves present (confirms fused views exist in this Variables)
nviews = Ref(0)
walk(x::NamedTuple) = foreach(k -> walk(getfield(x, k)), keys(x))
walk(x) = (x isa SpeedyWeather.LowerTriangularArrays.LowerTriangularArray || x isa SpeedyWeather.RingGrids.AbstractField) &&
    SpeedyWeather.is_view_entry(x) && (nviews[] += 1)
for g in SpeedyWeather.ALL_VARIABLE_GROUPS
    walk(getfield(vars, g))
end
println("view-backed leaves in Variables: ", nviews[])

@testset "to_vec(::Variables)" begin
    v, back = to_vec(vars)
    println("to_vec length = ", length(v), " eltype ", eltype(v))

    # round trip: reconstructing from the same vector reproduces the prognostic values
    rec = back(copy(v))
    @test rec isa SpeedyWeather.Variables
    @test to_vec(rec)[1] ≈ v

    # perturbation propagates: scaling the vector scales the reconstructed leaves (via the parents),
    # and the view leaves see it too (aliasing preserved by deepcopy)
    rec2 = back(2 .* v)
    @test to_vec(rec2)[1] ≈ 2 .* v
    # a prognostic view leaf in rec2 must reflect the 2x (i.e. it aliases the filled parent)
    vor = rec2.prognostic.vorticity
    vor0 = rec.prognostic.vorticity
    @test SpeedyWeather.is_view_entry(vor)            # it IS a view
    @test maximum(abs, vor.data) > 0
    @test vor.data ≈ 2 .* vor0.data                   # the view saw the parent update

    # full FD path on a cheap function returning a Variables — must not error (was Array->SubArray)
    f(x) = (y = deepcopy(x); y.prognostic.vorticity .*= 3; y)
    seed = deepcopy(vars)
    g = j′vp(central_fdm(5, 1), f, seed, vars)[1]
    @test g isa SpeedyWeather.Variables
    @test all(isfinite, to_vec(g)[1])
end
