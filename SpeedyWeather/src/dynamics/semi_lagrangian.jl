"""Semi-Lagrangian transport for the barotropic vorticity equation.

Along a trajectory the barotropic equation is `D(ζ+f)/Dt = sources`, so absolute vorticity is
materially conserved up to sources. Transport is then "evaluate the old field at the departure
point" rather than "evaluate a flux divergence", which removes the advective CFL restriction.

This file holds the geometry-facing half of the scheme: finding the departure points by backward
trajectory integration, and shifting absolute vorticity onto them via `RingGrids`' interpolation.
The time stepper that drives it lives in `time_stepping/steppers/semi_lagrangian.jl`."""

"""$(TYPEDSIGNATURES)
Convert a velocity in [m/s] and a time step in [s] into an angular displacement in [˚] on a
planet of `radius` [m]: `u*Δt` is a distance in [m], `/radius` makes it an angle in [rad],
`*180/π` converts to degrees. Stored on the time stepper so kernels take a single scalar."""
trajectory_time_step(Δt, radius) = Δt / radius * (180 / π)

"""$(TYPEDSIGNATURES)
Displace a point at `(lond, latd)` [˚E, ˚N] by velocities `u, v` [m/s] over the scaled time step
`Δt_deg` [˚*s/m] from [`trajectory_time_step`](@ref). Goes through `Particle` to reuse its tested
pole-crossing logic in `move` (crossing a pole flips the longitude by 180˚) and the `mod` wrapping
back into [0, 360˚E), [-90, 90˚N]."""
@inline function displace(lond::NF, latd::NF, u, v, Δt_deg) where {NF}
    dlat = v * Δt_deg

    # TODO: `cos(deg2rad(...))` rather than `cosd` to match `advect_2D`, see JuliaGPU/AMDGPU.jl#1041
    coslat = max(cos(deg2rad(latd)), eps(NF))   # prevents division by zero at the poles
    dlon = u * Δt_deg / coslat

    particle = mod(move(Particle{NF}(lond, latd), dlon, dlat))
    return particle.lon, particle.lat
end

"""$(TYPEDSIGNATURES)
Find the departure points of the backward trajectories arriving at every grid point.

Solves `x_d = x_a - Δt/2 (u*(x_a) + u*(x_d))` (the trapezoidal, or "iterated backward",
trajectory) by fixed-point iteration, starting from `x_d = x_a`. `u*` is the wind extrapolated to
the midpoint in time, see [`extrapolate_winds!`](@ref). Each iteration relocates the interpolation
stencil (`update_locator!`) and interpolates `u*, v*` onto the current estimate of the departure
points, so the cost is `n_iterations` locator updates plus `2*n_iterations` interpolations."""
function departure_points!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; locator, geometry) = time_stepping.interpolator
    (; n_iterations) = time_stepping
    Δt_deg = time_stepping.Δt_trajectory[]

    sl = vars.scratch.semi_lagrangian
    (; departure_lond, departure_latd, departure_u, departure_v) = sl

    # the (possibly time-extrapolated) wind that defines the trajectories
    u_star = sl.u_star
    v_star = sl.v_star
    extrapolate_winds!(vars, time_stepping, model)

    # arrival points are the grid points themselves
    arrival_lond, arrival_latd = time_stepping.arrival_lond, time_stepping.arrival_latd
    arch = architecture(departure_lond)
    npoints = length(departure_lond)

    # first guess: departure point = arrival point
    copyto!(departure_lond, arrival_lond)
    copyto!(departure_latd, arrival_latd)

    for _ in 1:n_iterations
        # relocate the interpolation stencil onto the current departure point estimate
        RingGrids.update_locator!(locator, geometry, departure_lond, departure_latd)

        # output is a plain vector, hitting `interpolate!(::AbstractVector, ::Field2D, ...)`
        RingGrids.interpolate!(departure_u, u_star, locator, geometry)
        RingGrids.interpolate!(departure_v, v_star, locator, geometry)

        launch!(
            arch, LinearWorkOrder, (npoints,), _departure_point_kernel!,
            departure_lond, departure_latd, arrival_lond, arrival_latd,
            u_star.data, v_star.data, departure_u, departure_v, Δt_deg
        )
    end
    return nothing
end

@kernel inbounds = true function _departure_point_kernel!(
        departure_lond, departure_latd, @Const(arrival_lond), @Const(arrival_latd),
        @Const(u_arrival), @Const(v_arrival), @Const(departure_u), @Const(departure_v), Δt_deg
    )
    ij = @index(Global, Linear)

    # trapezoidal rule: mean of the wind at the arrival and (current estimate of the) departure
    # point. Backward in time, hence the minus sign folded into -Δt_deg/2.
    u_mean = (u_arrival[ij] + departure_u[ij]) / 2
    v_mean = (v_arrival[ij] + departure_v[ij]) / 2

    lond, latd = displace(arrival_lond[ij], arrival_latd[ij], u_mean, v_mean, -Δt_deg / 2)
    departure_lond[ij] = lond
    departure_latd[ij] = latd
end

"""$(TYPEDSIGNATURES)
Fill `u_star, v_star` with the wind used to define the trajectories, extrapolated to the middle of
the time step. With `extrapolate_winds = true` this is the SETTLS-style `3/2 uⁿ - 1/2 uⁿ⁻¹`, second
order in time for a two-time-level scheme; otherwise just `uⁿ`, which is first order but avoids the
weakly-stable extrapolation. Falls back to `uⁿ` while `uⁿ⁻¹` is not yet available (first step)."""
function extrapolate_winds!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    sl = vars.scratch.semi_lagrangian
    (; clock) = vars.prognostic

    # grid u, v carry 2 steps for this time stepper: 1 = previous, 2 = current
    u_new = field_view(vars.grid.u, :, 1, 2)
    v_new = field_view(vars.grid.v, :, 1, 2)
    u_old = field_view(vars.grid.u, :, 1, 1)
    v_old = field_view(vars.grid.v, :, 1, 1)

    # no previous step available on the very first time step
    extrapolate = time_stepping.extrapolate_winds && clock.step_counter > 0
    w = extrapolate ? 3 // 2 : 1
    w_old = extrapolate ? -1 // 2 : 0

    @. sl.u_star = w * u_new + w_old * u_old
    @. sl.v_star = w * v_new + w_old * v_old
    return nothing
end

"""$(TYPEDSIGNATURES)
Shift absolute vorticity `ζ + f` to the departure points and subtract `f` again at the arrival
point, giving the transported relative vorticity `ζ*` on the grid. `f` is scaled by the same
`scale` as the (radius-scaled) prognostic vorticity, exactly as `_vorticity_flux_kernel!` does.

`departure_points!` must have run first: the locator still holds the departure-point stencil from
its last iteration, so no further `update_locator!` is needed here."""
function semi_lagrangian_transport!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; locator, geometry) = time_stepping.interpolator
    (; f) = model.coriolis
    scale = vars.prognostic.scale[]

    sl = vars.scratch.semi_lagrangian
    absolute_vorticity = sl.absolute_vorticity
    vorticity_departure = sl.vorticity_departure

    vor = field_view(vars.grid.vorticity, :, 1, 2)      # current grid vorticity
    (; whichring) = vor.grid
    arch = architecture(vor)

    # ζ + f on the grid, with f scaled on the fly as vorticity is
    launch!(
        arch, LinearWorkOrder, size(absolute_vorticity), _add_coriolis_kernel!,
        absolute_vorticity.data, vor.data, f, scale, whichring, true
    )

    # Shift to the departure points. NOTE: `.data` (a vector) and not the Field itself, because
    # `interpolate!(::Field, ::Field2D, ...)` short-circuits to a plain `copyto!` when the two
    # fields share a grid — which they do here, the departure points just are not the grid points.
    RingGrids.interpolate!(vorticity_departure.data, absolute_vorticity, locator, geometry)

    # and subtract f again at the arrival point
    launch!(
        arch, LinearWorkOrder, size(vorticity_departure), _add_coriolis_kernel!,
        vorticity_departure.data, vorticity_departure.data, f, scale, whichring, false
    )
    return nothing
end

@kernel inbounds = true function _add_coriolis_kernel!(
        out, @Const(vor), @Const(f), scale, @Const(whichring), add::Bool
    )
    ij = @index(Global, Linear)
    j = whichring[ij]
    out[ij] = vor[ij] + ifelse(add, f[j] * scale, -f[j] * scale)
end
