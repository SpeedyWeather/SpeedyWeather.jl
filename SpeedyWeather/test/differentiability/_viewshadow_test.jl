# Does a VIEW-PRESERVING zeroed shadow (deepcopy + zero) fix update_prognostic! on 1.12?
# vs the default make_zero (materializes views to Array -> Duplicated(view,Array) mismatch).
using SpeedyWeather, Enzyme, Test
import SpeedyWeather: is_view_entry, ALL_VARIABLE_GROUPS, FusedParent
import SpeedyWeather.LowerTriangularArrays: LowerTriangularArray
import SpeedyWeather.RingGrids: AbstractField
println("Julia ", VERSION)

# view-preserving zeroed shadow
zero_entry!(x::NamedTuple) = foreach(k -> zero_entry!(getfield(x,k)), keys(x))
zero_entry!(x::FusedParent) = zero_entry!(x.parent)
zero_entry!(x::LowerTriangularArray) = is_view_entry(x) || fill!(x.data, 0)
zero_entry!(x::AbstractField) = is_view_entry(x) || fill!(x.data, 0)
zero_entry!(::SubArray) = nothing
zero_entry!(x::AbstractArray) = eltype(x) <: Union{AbstractFloat,Complex} && fill!(x, 0)
zero_entry!(x) = nothing
function view_shadow(vars)
    z = deepcopy(vars)
    foreach(g -> zero_entry!(getfield(z, g)), ALL_VARIABLE_GROUPS)
    return z
end

spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
sim = initialize!(model); initialize!(sim)
vars = sim.variables
dvars = view_shadow(vars)
println("prognostic vorticity shadow type: ", typeof(dvars.prognostic.vorticity.data), " (view? ", is_view_entry(dvars.prognostic.vorticity), ")")
println("matches primal type? ", typeof(dvars.prognostic.vorticity) == typeof(vars.prognostic.vorticity))
dvars.prognostic.vorticity .= 1 + im
println(">>> entering autodiff update_prognostic! with VIEW-preserving shadow"); flush(stdout)
autodiff(set_runtime_activity(Reverse), SpeedyWeather.update_prognostic!, Const, Duplicated(vars, dvars), Const(model))
println(">>> autodiff RETURNED; nonzero: ", sum(to_vec(dvars)[1]) != 0)
