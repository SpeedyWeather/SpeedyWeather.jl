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
Velocities `u, v` [m/s] at `(lond, latd)` [˚E, ˚N] as a 3D Cartesian vector, tangent to the unit
sphere by construction. `u, v` are components in the *local* east/north frame, which rotates from
point to point, so two velocities at different points can only be combined (averaged, differenced)
after mapping them into this common frame — see [`departure_points!`](@ref)."""
@inline function cartesian_velocity(lond, latd, u, v)
    sinλ, cosλ = sincos(deg2rad(lond))
    sinφ, cosφ = sincos(deg2rad(latd))

    # east = (-sinλ, cosλ, 0), north = (-sinφcosλ, -sinφsinλ, cosφ)
    Vx = -u * sinλ - v * sinφ * cosλ
    Vy = u * cosλ - v * sinφ * sinλ
    Vz = v * cosφ
    return Vx, Vy, Vz
end

"""$(TYPEDSIGNATURES)
Displace a point at `(lond, latd)` [˚E, ˚N] by the 3D Cartesian velocity `V` [m/s] over the scaled
time step `Δt_over_radius` [s/m] from [`trajectory_time_step`](@ref), along a great circle.

Done in 3D Cartesian coordinates on the unit sphere rather than by incrementing longitude and
latitude. The lat/lon form needs `dlon = u*Δt/(radius*cos(lat))`, which is singular at the poles:
on an octahedral Gaussian T128 grid the outermost ring sits at 89.28˚N where `cos(lat) = 0.0125`,
so a 50 m/s wind over a 45 min step gives a 97˚ longitude jump — a straight line in (lon, lat) that
bears no relation to the actual trajectory. The Cartesian form has no coordinate singularity, is
exact for solid-body rotation, and follows the great circle rather than a rhumb-like path.

Backward trajectories are obtained by passing a negative `Δt_over_radius`."""
@inline function displace(lond::NF, latd::NF, Vx, Vy, Vz, Δt_over_radius) where {NF}
    sinλ, cosλ = sincos(deg2rad(lond))
    sinφ, cosφ = sincos(deg2rad(latd))
    rx, ry, rz = cosφ * cosλ, cosφ * sinλ, sinφ     # position on the unit sphere
    return rotate_along_great_circle(rx, ry, rz, Vx, Vy, Vz, Δt_over_radius)
end

"""$(TYPEDSIGNATURES)
Rotate the unit-sphere position `r` by the angle `|V| Δt_over_radius` along the great circle
spanned by `r` and `V`, returning `(lond, latd)` [˚E, ˚N]. The velocity does not have to be
tangent at `r`: its radial component is removed first, which matters because a mean of velocities
tangent at two *different* points is not itself tangent at either, and the rotation is only a
rotation for `V ⟂ r`."""
@inline function rotate_along_great_circle(rx::NF, ry, rz, Vx, Vy, Vz, Δt_over_radius) where {NF}
    Vr = Vx * rx + Vy * ry + Vz * rz
    Vx -= Vr * rx
    Vy -= Vr * ry
    Vz -= Vr * rz

    speed = sqrt(Vx^2 + Vy^2 + Vz^2)
    α = speed * Δt_over_radius                  # great-circle angle swept [rad], signed

    # rotate r by α along the great circle spanned by r and V:
    #   r_new = r*cos(α) + V̂*sin(α),  V̂ = V/speed
    # `sin(α)/speed` is written out so the speed → 0 limit (= Δt_over_radius) stays finite
    scale = ifelse(speed > eps(NF), sin(α) / speed, Δt_over_radius)
    cosα = cos(α)
    dx = rx * cosα + Vx * scale
    dy = ry * cosα + Vy * scale
    dz = rz * cosα + Vz * scale

    # back to degrees; asin puts latitude in [-90, 90] and the mod longitude in [0, 360)
    latd_new = rad2deg(asin(clamp(dz, -one(NF), one(NF))))
    lond_new = mod(rad2deg(atan(dy, dx)), 360)
    return lond_new, latd_new
end

"""$(TYPEDSIGNATURES)
One trapezoidal backward-trajectory iterate: the departure point reached from the arrival point
`(lond_a, latd_a)` [˚E, ˚N] by the mean of the arrival wind `u_a, v_a` [m/s], in local east/north
components, and the departure wind `Vd` [m/s], already in the common 3D frame.

Fuses [`cartesian_velocity`](@ref) and [`displace`](@ref) so that the arrival point's `sincos`
pair is evaluated once rather than once in each — this runs per grid point per iteration and the
transcendentals dominate it."""
@inline function departure_point(lond_a::NF, latd_a::NF, u_a, v_a, Vdx, Vdy, Vdz, Δt_over_radius) where {NF}
    sinλ, cosλ = sincos(deg2rad(lond_a))
    sinφ, cosφ = sincos(deg2rad(latd_a))
    rx, ry, rz = cosφ * cosλ, cosφ * sinλ, sinφ

    # the arrival wind in the common frame, reusing the sincos above (cf. `cartesian_velocity`)
    Vax = -u_a * sinλ - v_a * sinφ * cosλ
    Vay = u_a * cosλ - v_a * sinφ * sinλ
    Vaz = v_a * cosφ

    return rotate_along_great_circle(
        rx, ry, rz,
        (Vax + Vdx) / 2, (Vay + Vdy) / 2, (Vaz + Vdz) / 2,
        Δt_over_radius,
    )
end

"""$(TYPEDSIGNATURES)
Displace a point by the local-frame velocities `u, v` [m/s], see [`cartesian_velocity`](@ref)."""
@inline function displace(lond::NF, latd::NF, u, v, Δt_over_radius) where {NF}
    Vx, Vy, Vz = cartesian_velocity(lond, latd, u, v)
    return displace(lond, latd, Vx, Vy, Vz, Δt_over_radius)
end

"""$(TYPEDSIGNATURES)
Find the departure points of the backward trajectories arriving at every grid point.

Solves `x_d = x_a - Δt/2 (u*(x_a) + u*(x_d))` (the trapezoidal, or "iterated backward",
trajectory) by fixed-point iteration, starting from `x_d = x_a`. `u*` is the wind extrapolated to
the midpoint in time, see [`extrapolate_winds!`](@ref). Each iteration relocates the interpolation
stencil (`update_locator!`) and interpolates `u*, v*` onto the current estimate of the departure
points. The two winds entering the mean are mapped to a common 3D Cartesian frame first, see
[`cartesian_velocity`](@ref): `u, v` are components in the local east/north frame, and adding
them across two different points is otherwise wrong by the meridian convergence — which grows
towards the poles until the iteration diverges rather than converges.

The first guess needs no interpolation at all (the departure point starts at the arrival point,
where the wind is already known) and the winds at the final departure points are never read, so the
cost is `n_iterations - 1` locator updates plus `2*(n_iterations - 1)` interpolations — for the
default `n_iterations = 2` that is one of each pair, not two.

These use `trajectory_locator`, which only has to be accurate enough to place the departure point;
the damping that matters accumulates on the transported field, interpolated once per step by
`semi_lagrangian_transport!` with the more expensive `locator`.

The arrival points are the grid points, taken straight from `model.geometry`."""
function departure_points!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; geometry, n_iterations) = time_stepping
    Δt_over_radius = time_stepping.Δt_trajectory[]

    sl = vars.dynamics.semi_lagrangian
    (; departure_lond, departure_latd, departure_u, departure_v, u_trajectory, v_trajectory) = sl
    locator = sl.trajectory_locator

    # the (possibly time-extrapolated) wind that defines the trajectories
    extrapolate_winds!(vars, time_stepping, model)

    # arrival points are the grid points themselves
    arrival_lond, arrival_latd = model.geometry.londs, model.geometry.latds
    arch = architecture(departure_lond)
    npoints = length(departure_lond)

    # First guess for the departure point is the arrival point, where the wind is just `u_trajectory` at
    # the grid point: interpolating there would be the identity. So seed the departure point and
    # the winds there directly rather than paying a locator update and two interpolations.
    copyto!(departure_lond.data, arrival_lond)
    copyto!(departure_latd.data, arrival_latd)
    copyto!(departure_u.data, u_trajectory.data)
    copyto!(departure_v.data, v_trajectory.data)

    for iteration in 1:n_iterations
        launch!(
            arch, LinearWorkOrder, (npoints,), _departure_point_kernel!,
            departure_lond.data, departure_latd.data, arrival_lond, arrival_latd,
            u_trajectory.data, v_trajectory.data, departure_u.data, departure_v.data, Δt_over_radius
        )

        # the winds at the *final* departure points are never read — the next thing that happens
        # is the tracer shift, which needs the points, not the winds — so skip the last refresh
        iteration == n_iterations && break

        # relocate the trajectory stencil onto the current departure point estimate.
        # NOTE: `.data` (a vector) throughout, because `interpolate!(::Field, ::Field2D, ...)`
        # short-circuits to a plain `copyto!` when the two fields share a grid — which they do,
        # the departure points just are not the grid points.
        RingGrids.update_locator!(locator, geometry, departure_lond.data, departure_latd.data)
        RingGrids.interpolate!(departure_u.data, u_trajectory, locator, geometry)
        RingGrids.interpolate!(departure_v.data, v_trajectory, locator, geometry)
    end
    return nothing
end

@kernel inbounds = true function _departure_point_kernel!(
        departure_lond, departure_latd, @Const(arrival_lond), @Const(arrival_latd),
        @Const(u_arrival), @Const(v_arrival), @Const(departure_u), @Const(departure_v), Δt_over_radius
    )
    ij = @index(Global, Linear)

    # trapezoidal rule: mean of the wind at the arrival and (current estimate of the) departure
    # point. The two live in different local east/north frames, so they are mapped to a common
    # 3D Cartesian frame before being averaged — adding the components directly is only valid if
    # the frames coincide, and the mismatch is the meridian convergence, O(Δλ sinφ), which is
    # small in midlatitudes but O(1) next to the poles where a trajectory spans a large Δλ.
    Vd = cartesian_velocity(departure_lond[ij], departure_latd[ij], departure_u[ij], departure_v[ij])

    # the mean wind is already the 1/2 in x_d = x_a - Δt/2 (V(x_a) + V(x_d)), so the full time
    # step is used here; backward in time, hence the minus sign
    lond, latd = departure_point(
        arrival_lond[ij], arrival_latd[ij], u_arrival[ij], v_arrival[ij],
        Vd..., -Δt_over_radius
    )
    departure_lond[ij] = lond
    departure_latd[ij] = latd
end

"""$(TYPEDSIGNATURES)
Fill `u_trajectory, v_trajectory` with the wind used to define the trajectories, extrapolated to
the middle of the time step. With `extrapolate_winds = true` this is the SETTLS-style
`3/2 uⁿ - 1/2 uⁿ⁻¹`, second order in time for a two-time-level scheme; otherwise just `uⁿ`, which
is first order.

This is not a free choice: the trajectory carries the Rossby term via `f(x_d) - f(x_a)`, so a wind
at `tⁿ` makes that term forward Euler and hence unconditionally unstable. `extrapolate_winds =
false` is a diagnostic only, it blows up within days at a 1 h time step.

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

    # Linear extrapolation in time from uⁿ⁻¹, uⁿ to the middle of the step, which is where the
    # trapezoidal trajectory wants its wind:
    #   u(tⁿ + Δt/2) ≈ uⁿ + (Δt/2) (uⁿ - uⁿ⁻¹)/Δt = 3/2 uⁿ - 1/2 uⁿ⁻¹
    # hence the weights. Without it the wind stays at tⁿ, one step behind (see the docstring:
    # that is forward Euler on the Rossby term and unstable, not merely first order).
    NF = eltype(sl.u_trajectory)
    w, w_old = time_stepping.extrapolate_winds ? (NF(3 // 2), NF(-1 // 2)) : (one(NF), zero(NF))

    @. sl.u_trajectory = w * u_new + w_old * u_old
    @. sl.v_trajectory = w * v_new + w_old * v_old
    return nothing
end

"""$(TYPEDSIGNATURES)
Shift absolute vorticity `ζ + f` to the departure points and subtract `f` again at the arrival
point, giving the transported relative vorticity `ζ*` on the grid. `f` is scaled by the same
`scale` as the (radius-scaled) prognostic vorticity, exactly as `_vorticity_flux_kernel!` does.

`ζ + f` is formed in place in the grid vorticity rather than in a dedicated array: nothing reads
`vars.grid.vorticity` again before the closing `transform!` of the time step overwrites it from
the new spectral state.

`departure_points!` must have run first to fill the departure points; the field locator's own
stencil is located here, as it is generally a different (higher order) one than the trajectory's."""
function semi_lagrangian_transport!(
        vars::Variables,
        time_stepping::AbstractSemiLagrangian,
        model::AbstractModel,
    )
    (; geometry) = time_stepping
    (; f) = model.coriolis
    scale = vars.prognostic.scale[]

    sl = vars.dynamics.semi_lagrangian
    vorticity_departure = sl.vorticity_departure
    locator = sl.locator

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
