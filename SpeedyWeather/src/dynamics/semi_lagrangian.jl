"""Semi-Lagrangian transport for the barotropic vorticity equation.

Along a trajectory the barotropic equation is `D(ζ+f)/Dt = sources`, so absolute vorticity is
materially conserved up to sources. Transport is then "evaluate the old field at the departure
point" rather than "evaluate a flux divergence", which removes the advective CFL restriction.

This file holds the geometry-facing half of the scheme: finding the departure points by backward
trajectory integration, and shifting absolute vorticity onto them via `RingGrids`' interpolation.
The time stepper that drives it lives in `time_stepping/steppers/semi_lagrangian.jl`."""

"""$(TYPEDSIGNATURES)
Scale a time step [s] for use in [`displace`](@ref): `u*Δt/radius` is the great-circle angle
[rad] a point at velocity `u` [m/s] sweeps in time `Δt` on a planet of `radius` [m]."""
trajectory_time_step(Δt, radius) = Δt / radius

"""$(TYPEDSIGNATURES)
Displace a point at `(lond, latd)` [˚E, ˚N] by velocities `u, v` [m/s] over the scaled time step
`Δt_over_radius` [s/m] from [`trajectory_time_step`](@ref), along a great circle.

Done in 3D Cartesian coordinates on the unit sphere rather than by incrementing longitude and
latitude. The lat/lon form needs `dlon = u*Δt/(radius*cos(lat))`, which is singular at the poles:
on an octahedral Gaussian T128 grid the outermost ring sits at 89.28˚N where `cos(lat) = 0.0125`,
so a 50 m/s wind over a 45 min step gives a 97˚ longitude jump — a straight line in (lon, lat) that
bears no relation to the actual trajectory. The Cartesian form has no coordinate singularity, is
exact for solid-body rotation, and follows the great circle rather than a rhumb-like path.

Backward trajectories are obtained by passing a negative `Δt_over_radius`."""
@inline function displace(lond::NF, latd::NF, u, v, Δt_over_radius) where {NF}
    sinλ, cosλ = sincos(deg2rad(lond))
    sinφ, cosφ = sincos(deg2rad(latd))

    # position on the unit sphere and the local east/north unit vectors there
    rx, ry, rz = cosφ * cosλ, cosφ * sinλ, sinφ
    ex, ey = -sinλ, cosλ                        # east, ez = 0
    nx, ny, nz = -sinφ * cosλ, -sinφ * sinλ, cosφ   # north

    # velocity as a 3D vector, tangent to the sphere by construction
    Vx = u * ex + v * nx
    Vy = u * ey + v * ny
    Vz = v * nz                                 # u * ez = 0

    speed = sqrt(Vx^2 + Vy^2 + Vz^2)
    α = speed * Δt_over_radius                  # great-circle angle swept [rad], signed

    # rotate r by α along the great circle spanned by r and V:
    #   r_new = r*cos(α) + V̂*sin(α),  V̂ = V/speed
    # `sin(α)/speed` is written out so the speed → 0 limit (= Δt_over_radius) stays finite
    scale = ifelse(speed > eps(NF), sin(α) / speed, Δt_over_radius)
    dx = rx * cos(α) + Vx * scale
    dy = ry * cos(α) + Vy * scale
    dz = rz * cos(α) + Vz * scale

    # back to degrees; asin puts latitude in [-90, 90] and the mod longitude in [0, 360)
    latd_new = rad2deg(asin(clamp(dz, -one(NF), one(NF))))
    lond_new = mod(rad2deg(atan(dy, dx)), 360)
    return lond_new, latd_new
end

"""$(TYPEDSIGNATURES)
Find the departure points of the backward trajectories arriving at every grid point.

Solves `x_d = x_a - Δt/2 (u*(x_a) + u*(x_d))` (the trapezoidal, or "iterated backward",
trajectory) by fixed-point iteration, starting from `x_d = x_a`. `u*` is the wind extrapolated to
the midpoint in time, see [`extrapolate_winds!`](@ref). Each iteration relocates the interpolation
stencil (`update_locator!`) and interpolates `u*, v*` onto the current estimate of the departure
points. The first guess needs no interpolation at all (the departure point starts at the arrival
point, where the wind is already known) and the winds at the final departure points are never read,
so the cost is `n_iterations - 1` locator updates plus `2*(n_iterations - 1)` interpolations — for
the default `n_iterations = 2` that is one of each pair, not two.

These use `trajectory_interpolator`, which only has to be accurate enough to place the departure
point; the damping that matters accumulates on the transported field, interpolated once per step by
`semi_lagrangian_transport!` with the more expensive `interpolator`.

The arrival points are the grid points, taken straight from `model.geometry`."""
function departure_points!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; locator, geometry) = time_stepping.trajectory_interpolator
    (; n_iterations) = time_stepping
    Δt_over_radius = time_stepping.Δt_trajectory[]

    sl = vars.dynamics.semi_lagrangian
    (; departure_lond, departure_latd, departure_u, departure_v, u_star, v_star) = sl

    # the (possibly time-extrapolated) wind that defines the trajectories
    extrapolate_winds!(vars, time_stepping, model)

    # arrival points are the grid points themselves
    arrival_lond, arrival_latd = model.geometry.londs, model.geometry.latds
    arch = architecture(departure_lond)
    npoints = length(departure_lond)

    # First guess for the departure point is the arrival point, where the wind is just `u_star` at
    # the grid point: interpolating there would be the identity. So seed the departure winds
    # directly rather than paying a locator update and two interpolations to recompute them.
    copyto!(departure_u.data, u_star.data)
    copyto!(departure_v.data, v_star.data)

    for iteration in 1:n_iterations
        launch!(
            arch, LinearWorkOrder, (npoints,), _departure_point_kernel!,
            departure_lond.data, departure_latd.data, arrival_lond, arrival_latd,
            u_star.data, v_star.data, departure_u.data, departure_v.data, Δt_over_radius
        )

        # the winds at the *final* departure points are never read — the next thing that happens
        # is the tracer shift, which needs the points, not the winds — so skip the last refresh
        iteration == n_iterations && break

        # relocate the trajectory stencil onto the current departure point estimate.
        # NOTE: `.data` (a vector) throughout, because `interpolate!(::Field, ::Field2D, ...)`
        # short-circuits to a plain `copyto!` when the two fields share a grid — which they do,
        # the departure points just are not the grid points.
        RingGrids.update_locator!(locator, geometry, departure_lond.data, departure_latd.data)
        RingGrids.interpolate!(departure_u.data, u_star, locator, geometry)
        RingGrids.interpolate!(departure_v.data, v_star, locator, geometry)
    end
    return nothing
end

@kernel inbounds = true function _departure_point_kernel!(
        departure_lond, departure_latd, @Const(arrival_lond), @Const(arrival_latd),
        @Const(u_arrival), @Const(v_arrival), @Const(departure_u), @Const(departure_v), Δt_over_radius
    )
    ij = @index(Global, Linear)

    # trapezoidal rule: mean of the wind at the arrival and (current estimate of the) departure
    # point. Backward in time, hence the minus sign folded into -Δt_over_radius/2.
    u_mean = (u_arrival[ij] + departure_u[ij]) / 2
    v_mean = (v_arrival[ij] + departure_v[ij]) / 2

    lond, latd = displace(arrival_lond[ij], arrival_latd[ij], u_mean, v_mean, -Δt_over_radius / 2)
    departure_lond[ij] = lond
    departure_latd[ij] = latd
end

"""$(TYPEDSIGNATURES)
Fill `u_star, v_star` with the wind used to define the trajectories, extrapolated to the middle of
the time step. With `extrapolate_winds = true` this is the SETTLS-style `3/2 uⁿ - 1/2 uⁿ⁻¹`, second
order in time for a two-time-level scheme; otherwise just `uⁿ`, which is first order but avoids the
weakly-stable extrapolation.

No special case for the first time step: `move_prognostic_grid_variables_back!` is called by the
barotropic `transform!` with `initialize=true`, so `uⁿ⁻¹ = uⁿ` going into the first step and the
extrapolation reduces to `uⁿ` on its own."""
function extrapolate_winds!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    sl = vars.dynamics.semi_lagrangian

    # grid u, v carry 2 steps for this time stepper: 1 = previous, 2 = current
    u_new = field_view(vars.grid.u, :, 1, 2)
    v_new = field_view(vars.grid.v, :, 1, 2)
    u_old = field_view(vars.grid.u, :, 1, 1)
    v_old = field_view(vars.grid.v, :, 1, 1)

    NF = eltype(sl.u_star)
    w, w_old = time_stepping.extrapolate_winds ? (NF(3 // 2), NF(-1 // 2)) : (one(NF), zero(NF))

    @. sl.u_star = w * u_new + w_old * u_old
    @. sl.v_star = w * v_new + w_old * v_old
    return nothing
end

"""$(TYPEDSIGNATURES)
Shift absolute vorticity `ζ + f` to the departure points and subtract `f` again at the arrival
point, giving the transported relative vorticity `ζ*` on the grid. `f` is scaled by the same
`scale` as the (radius-scaled) prognostic vorticity, exactly as `_vorticity_flux_kernel!` does.

`ζ + f` is formed in place in the grid vorticity rather than in a dedicated array: nothing reads
`vars.grid.vorticity` again before the closing `transform!` of the time step overwrites it from
the new spectral state.

`departure_points!` must have run first to fill the departure points; the field interpolator's own
stencil is located here, as it is generally a different (higher order) one than the trajectory's."""
function semi_lagrangian_transport!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; locator, geometry) = time_stepping.interpolator
    (; f) = model.coriolis
    scale = vars.prognostic.scale[]

    sl = vars.dynamics.semi_lagrangian
    vorticity_departure = sl.vorticity_departure

    # locate the field stencil on the departure points found by `departure_points!`
    RingGrids.update_locator!(locator, geometry, sl.departure_lond.data, sl.departure_latd.data)
    vor = field_view(vars.grid.vorticity, :, 1, 2)      # current grid vorticity
    (; whichring) = vor.grid
    arch = architecture(vor)

    # ζ += f in place, f scaled on the fly as vorticity is
    launch!(
        arch, LinearWorkOrder, size(vor), _add_coriolis_kernel!,
        vor.data, f, scale, whichring, true
    )

    # shift to the departure points, `.data` to hit the vector method (see departure_points!)
    RingGrids.interpolate!(vorticity_departure.data, vor, locator, geometry)

    # and subtract f again at the arrival point
    launch!(
        arch, LinearWorkOrder, size(vorticity_departure), _add_coriolis_kernel!,
        vorticity_departure.data, f, scale, whichring, false
    )
    return nothing
end

@kernel inbounds = true function _add_coriolis_kernel!(
        vor, @Const(f), scale, @Const(whichring), add::Bool
    )
    ij = @index(Global, Linear)
    j = whichring[ij]
    vor[ij] = vor[ij] + ifelse(add, f[j] * scale, -f[j] * scale)
end
