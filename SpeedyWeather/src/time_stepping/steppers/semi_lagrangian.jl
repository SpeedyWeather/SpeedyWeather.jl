export SemiLagrangian

"""Two-time-level semi-Lagrangian time stepping for the barotropic vorticity equation.

Transport is done by shifting absolute vorticity `ζ+f` to the departure points of the backward
trajectories (`dynamics/semi_lagrangian.jl`) rather than through the vorticity flux, removing the
advective CFL restriction. The remaining terms are integrated as an exponential integrator (ETD1):
sources explicitly, hyperdiffusion exactly through `φ₁`. One step is

    ζⁿ⁺¹ = exp(z) ζ* + Δt φ₁(z) S,    z = Δt∇²ⁿ,  ζ* = ζⁿ shifted to the departure points

with `S` the source tendency (curl of the forcing, drag). This is a Lie–Trotter split of two
operators each handled exactly; the splitting error is `O(Δt²)` as the two do not commute.

Only `BarotropicModel` is supported.
$(TYPEDFIELDS)"""
mutable struct SemiLagrangian{NF, S, B, MS, IP, VT} <: AbstractSemiLagrangian
    "[OPTION] Time step for T32, scale linearly to spectral resolution `truncation`"
    Δt_at_T32::S

    "[OPTION] Adjust `Δt_at_T32` with the output `interval` to output exactly after integer time steps"
    adjust_with_output::B

    "[OPTION] Fixed-point iterations for the backward trajectory, 2 is usually enough"
    n_iterations::Int

    "[OPTION] Extrapolate the trajectory wind in time as 3/2 uⁿ - 1/2 uⁿ⁻¹ (SETTLS-style)"
    extrapolate_winds::Bool

    "[OPTION] Use the exact exponential (ETD1) for hyperdiffusion, set on `model.horizontal_diffusion`"
    exponential_diffusion::Bool

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt::NF

    "[DERIVED] Δt scaled to convert [m/s] into [˚], see `trajectory_time_step`"
    Δt_trajectory::Base.RefValue{NF}

    "[DERIVED] Interpolator (grid geometry + locator) used to evaluate fields at departure points"
    interpolator::IP

    "[DERIVED] Longitudes [˚E] of the arrival points, i.e. of the grid points"
    arrival_lond::VT

    "[DERIVED] Latitudes [˚N] of the arrival points, i.e. of the grid points"
    arrival_latd::VT
end

"""$(TYPEDSIGNATURES)
Generator function for `SemiLagrangian` using `spectral_grid` for resolution."""
function SemiLagrangian(
        spectral_grid::SpectralGrid;
        Δt_at_T32 = Minute(60),         # larger than Leapfrog's default: no advective CFL limit
        adjust_with_output = true,
        n_iterations = 2,
        extrapolate_winds = true,
        exponential_diffusion = true,
    )
    (; NF, truncation, grid) = spectral_grid

    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T32), truncation, DEFAULT_RADIUS, adjust_with_output)
    Δt::NF = Δt_millisec.value / 1000

    # one interpolation target per grid point: the departure point of the trajectory arriving there
    npoints = RingGrids.get_npoints(grid)
    interpolator = RingGrids.AnvilInterpolator(grid, npoints; NF)

    londs, latds = RingGrids.get_londlatds(grid)
    arrival_lond = on_architecture(architecture(grid), NF.(londs))
    arrival_latd = on_architecture(architecture(grid), NF.(latds))

    return SemiLagrangian{NF, Second, Bool, Millisecond, typeof(interpolator), typeof(arrival_lond)}(
        Second(Δt_at_T32), adjust_with_output, n_iterations, extrapolate_winds,
        exponential_diffusion, Δt_millisec, Δt, Ref(zero(NF)),
        interpolator, arrival_lond, arrival_latd,
    )
end

function initialize!(L::SemiLagrangian, model::AbstractModel)
    model isa Barotropic || throw(ArgumentError(
            "SemiLagrangian time stepping is currently only implemented for BarotropicModel, got $(typeof(model))."
        ))

    calculate_Δt!(L, model)
    L.Δt_trajectory[] = trajectory_time_step(L.Δt, model.planet.radius)

    # The ETD1 update below is only exact for the diffusion if `impl` holds φ₁ rather than the
    # backward-Euler factor. `initialize!(model.time_stepping, ...)` runs before
    # `initialize!(model.horizontal_diffusion, ...)` in `initialize!(::Barotropic)`, so setting
    # the flag here still takes effect when the diffusion arrays are precomputed.
    if L.exponential_diffusion && hasproperty(model.horizontal_diffusion, :exponential)
        model.horizontal_diffusion.exponential = true
    end
    return nothing
end

# HOW MANY STEPS DO VARIABLES NEED?
# two-time-level: a single spectral state ...
prognostic_spectral_steps(::AbstractSemiLagrangian) = 1
# ... but two grid time levels of u, v for the SETTLS wind extrapolation
prognostic_grid_steps(::AbstractSemiLagrangian, ::Barotropic) = 2
tendency_steps(::AbstractSemiLagrangian) = 1

# WHICH STEP TO READ WHEN
# spectral variables have a single step
@inline which_prognostic_step(var, ::AbstractSemiLagrangian, ::STEP_COMPONENT) = 1
# grid variables have two, the 2nd is the current one (the 1st is the previous time step)
@inline which_prognostic_step(var::AbstractField, ::AbstractSemiLagrangian, ::STEP_COMPONENT) = 2

# the whole point: the trajectory shift replaces the vorticity flux
@inline advection_factor(::AbstractSemiLagrangian) = 0

# diffusion is applied to the transported state, not corrected implicitly alongside
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::Nothing, ::AbstractSemiLagrangian) = true

"""$(TYPEDSIGNATURES) Work arrays for the departure point search and the trajectory shift."""
function variables(L::SemiLagrangian, model::AbstractModel)
    npoints = RingGrids.get_npoints(model.spectral_grid.grid)
    ns = :semi_lagrangian
    return (
        ScratchVariable(:departure_lond, VectorDim(npoints), namespace = ns, desc = "Departure point longitude", units = "˚E"),
        ScratchVariable(:departure_latd, VectorDim(npoints), namespace = ns, desc = "Departure point latitude", units = "˚N"),
        ScratchVariable(:departure_u, VectorDim(npoints), namespace = ns, desc = "Zonal wind at departure point", units = "m/s"),
        ScratchVariable(:departure_v, VectorDim(npoints), namespace = ns, desc = "Meridional wind at departure point", units = "m/s"),
        ScratchVariable(:u_star, Grid2D(), namespace = ns, desc = "Time-extrapolated zonal wind for trajectories", units = "m/s"),
        ScratchVariable(:v_star, Grid2D(), namespace = ns, desc = "Time-extrapolated meridional wind for trajectories", units = "m/s"),
        ScratchVariable(:absolute_vorticity, Grid2D(), namespace = ns, desc = "Absolute vorticity ζ+f", units = "1/s"),
        ScratchVariable(:vorticity_departure, Grid2D(), namespace = ns, desc = "Vorticity shifted to departure points", units = "1/s"),
    )
end

"""$(TYPEDSIGNATURES)
Retain the current grid `u, v` (step 2) as the previous ones (step 1) for the next step's wind
extrapolation. The barotropic `transform!` does not call this hook (its default is a no-op and
`Leapfrog` uses a single grid step in 2D), so the semi-Lagrangian step calls it itself."""
function move_prognostic_grid_variables_back!(
        vars::Variables,
        ::AbstractSemiLagrangian,
        ::Barotropic,
    )
    copy_step_back!(parent(vars.fused.uv_grid))
    return nothing
end

"""$(TYPEDSIGNATURES)
One semi-Lagrangian time step for the barotropic model.

The ordering matters: `semi_lagrangian_transport!` writes the transported state `ζ*` into
`vars.prognostic.vorticity` *before* `horizontal_diffusion!` reads it, so the diffusion's
`(tendency + expl*var)*impl` form picks up the already-transported field and the closing Euler
update becomes exactly `ζⁿ⁺¹ = exp(z)ζ* + Δt φ₁(z) S`."""
function time_step!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::Barotropic,
    )
    (!isnothing(model.feedback) && model.feedback.nans_detected) && return nothing
    reset_tendencies!(vars, time_stepping)

    # SOURCES. `advection_factor(::AbstractSemiLagrangian) = 0` switches the vorticity flux off,
    # so this reduces to S = ∇×(Fᵤ, Fᵥ) - cζ in spectral space, radius-scaled as usual.
    dynamics_tendencies!(vars, model)

    # TRANSPORT. Backward trajectories, then shift ζ+f onto their departure points.
    departure_points!(vars, time_stepping, model)
    semi_lagrangian_transport!(vars, time_stepping, model)

    # retain the current grid u, v as the previous ones before `transform!` overwrites them
    move_prognostic_grid_variables_back!(vars, time_stepping, model)

    # ζ* back to spectral, overwriting the prognostic state
    vorticity = get_step(vars.prognostic.vorticity, 1)
    transform!(
        vorticity, vars.scratch.semi_lagrangian.vorticity_departure,
        vars.scratch.transform_memory, model.spectral_transform
    )

    # DIFFUSION + SOURCES, as an ETD1 step on the transported state
    horizontal_diffusion!(vars, model)
    update_prognostic!(vars, model)

    transform!(vars, model)
    particle_advection!(vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Forward Euler on the spectral state. Correct here because the transport is already baked into the
state by the trajectory shift and the diffusion into `φ₁` by `horizontal_diffusion!`."""
function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        clock::Clock,
        time_stepping::AbstractSemiLagrangian,
        implicit::Union{Nothing, AbstractImplicit},
        ::AbstractModel,
        scale::Real = 1,
    )
    (; Δt) = time_stepping
    Δt /= oftype(Δt, scale)

    var_step = get_step(var, 1)
    var_tend = get_tendency_step(tendency, time_stepping, time_stepping)

    launch!(
        architecture(var_tend), SpectralWorkOrder, size(var_tend), semi_lagrangian_kernel!,
        var_step, var_tend, Δt
    )
    return nothing
end

@kernel inbounds = true function semi_lagrangian_kernel!(var, @Const(tendency), Δt)
    lmk = @index(Global, Linear)
    var[lmk] = var[lmk] + Δt * tendency[lmk]
end
