export interpolate_3D!, Center, Face

"""Marker for fields located on full vertical levels (cell centers): u, v, T, etc."""
struct Center end

"""Marker for fields located on half vertical levels (lower cell face, k+½): vertical
velocity w. The top face (σ=0, k=½) is not stored; it is implicitly zero by the
kinematic boundary condition. The bottom face (σ=1) is explicitly stored as zero."""
struct Face end

# Zero-gradient (von Neumann) BC: σ outside [σ_levels[1], σ_levels[end]]
# returns the boundary value unchanged (gradient assumed zero at boundary).
@inline function find_vertical_bracket(σ, σ_levels)
    n = length(σ_levels)
    k_lo = 1
    for k in 1:(n - 1)
        σ_levels[k] <= σ && (k_lo = k)
    end
    k_hi = k_lo + 1
    Δσ = σ_levels[k_hi] - σ_levels[k_lo]
    α = clamp((σ - σ_levels[k_lo]) / Δσ, zero(σ), one(σ))
    return k_lo, k_hi, α
end

# Compute north and south pole ring averages per vertical layer on device.
@kernel inbounds = true function _compute_pole_averages_kernel!(
        north_pole_average, south_pole_average, A_data, ring_starts, nlons, nlat
    )
    k = @index(Global, Linear)
    n_north = nlons[1]
    rs_north = ring_starts[1]
    north_sum = zero(eltype(north_pole_average))
    for i in 0:(n_north - 1)
        north_sum += A_data[rs_north + i, k]
    end
    north_pole_average[k] = north_sum / n_north
    n_south = nlons[nlat]
    rs_south = ring_starts[nlat]
    south_sum = zero(eltype(south_pole_average))
    for i in 0:(n_south - 1)
        south_sum += A_data[rs_south + i, k]
    end
    south_pole_average[k] = south_sum / n_south
end

# Full-level (Center) interpolation kernel.
@kernel inbounds = true function _interpolate_3D_kernel!(
        Aout, A_data, locator, positions, σ_levels,
    )
    (; ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys, north_pole_average, south_pole_average) = locator
    i = @index(Global, Linear)
    k_lo, k_hi, α = find_vertical_bracket(positions[i].σ, σ_levels)

    a = ij_as[i] == 0 ? north_pole_average[k_lo] : A_data[ij_as[i], k_lo]
    b = ij_as[i] == 0 ? north_pole_average[k_lo] : A_data[ij_bs[i], k_lo]
    c = ij_cs[i] == -1 ? south_pole_average[k_lo] : A_data[ij_cs[i], k_lo]
    d = ij_cs[i] == -1 ? south_pole_average[k_lo] : A_data[ij_ds[i], k_lo]
    ab_lo = a + (b - a) * Δabs[i]
    val_lo = ab_lo + (c + (d - c) * Δcds[i] - ab_lo) * Δys[i]

    a = ij_as[i] == 0 ? north_pole_average[k_hi] : A_data[ij_as[i], k_hi]
    b = ij_as[i] == 0 ? north_pole_average[k_hi] : A_data[ij_bs[i], k_hi]
    c = ij_cs[i] == -1 ? south_pole_average[k_hi] : A_data[ij_cs[i], k_hi]
    d = ij_cs[i] == -1 ? south_pole_average[k_hi] : A_data[ij_ds[i], k_hi]
    ab_hi = a + (b - a) * Δabs[i]
    val_hi = ab_hi + (c + (d - c) * Δcds[i] - ab_hi) * Δys[i]

    Aout[i] = val_lo + (val_hi - val_lo) * α
end

# Half-level (Face) interpolation kernel for fields like vertical velocity w.
# Layer k in A_data stores the value at the lower face of layer k (i.e. σ at k+½).
# The top face (σ=0, k=½) is not stored; it is implicitly zero (kinematic BC).
# The bracket indices into σ_levels_half give k_lo, k_hi; the data index is k_lo-1
# (the face above k_lo). When k_data_lo == 0 the upper neighbour is the implicit
# top BC → value is zero.
@kernel inbounds = true function _interpolate_3D_face_kernel!(
        Aout, A_data, locator, positions, σ_levels,
    )
    (; ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys, north_pole_average, south_pole_average) = locator
    i = @index(Global, Linear)
    k_lo, k_hi, α = find_vertical_bracket(positions[i].σ, σ_levels)

    # k_lo - 1 is the data index for the upper bracket face;
    # 0 means the implicit top BC (σ=0, w=0), not stored in A_data.
    k_data_lo = k_lo - 1
    k_data_hi = k_hi - 1   # always >= 1 since k_hi = k_lo + 1 >= 2

    Z = zero(eltype(A_data))

    # val_lo: horizontal interpolation at the upper face of the bracket
    val_lo = if k_data_lo == 0
        Z
    else
        a = ij_as[i] == 0 ? north_pole_average[k_data_lo] : A_data[ij_as[i], k_data_lo]
        b = ij_as[i] == 0 ? north_pole_average[k_data_lo] : A_data[ij_bs[i], k_data_lo]
        c = ij_cs[i] == -1 ? south_pole_average[k_data_lo] : A_data[ij_cs[i], k_data_lo]
        d = ij_cs[i] == -1 ? south_pole_average[k_data_lo] : A_data[ij_ds[i], k_data_lo]
        ab = a + (b - a) * Δabs[i]
        ab + (c + (d - c) * Δcds[i] - ab) * Δys[i]
    end

    # val_hi: horizontal interpolation at the lower face of the bracket — always in storage
    a = ij_as[i] == 0 ? north_pole_average[k_data_hi] : A_data[ij_as[i], k_data_hi]
    b = ij_as[i] == 0 ? north_pole_average[k_data_hi] : A_data[ij_bs[i], k_data_hi]
    c = ij_cs[i] == -1 ? south_pole_average[k_data_hi] : A_data[ij_cs[i], k_data_hi]
    d = ij_cs[i] == -1 ? south_pole_average[k_data_hi] : A_data[ij_ds[i], k_data_hi]
    ab = a + (b - a) * Δabs[i]
    val_hi = ab + (c + (d - c) * Δcds[i] - ab) * Δys[i]

    Aout[i] = val_lo + (val_hi - val_lo) * α
end

"""$(TYPEDSIGNATURES)
Vertically blended anvil interpolation at each particle's σ coordinate for fields on
full vertical levels (cell centers), e.g. u, v, temperature. Pole ring averages are
stored in `locator.north_pole_average` and `locator.south_pole_average`
(pre-allocated via the 3D `AnvilLocator` constructor) and reused across calls
without CPU allocations. Pass `σ_levels_full` from `model.geometry`."""
function interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels, ::Center)
    (; npoints_output, north_pole_average, south_pole_average) = locator
    A_data = A.data
    arch = architecture(Aout)
    (; ring_starts, nlons, nlat) = geometry
    launch!(arch, LinearWorkOrder, (size(A_data, 2),), _compute_pole_averages_kernel!,
            north_pole_average, south_pole_average, A_data, ring_starts, nlons, nlat)
    launch!(
        arch, LinearWorkOrder, (npoints_output,),
        _interpolate_3D_kernel!,
        Aout, A_data, locator, positions, σ_levels,
    )
    return Aout
end

"""$(TYPEDSIGNATURES)
Vertically blended anvil interpolation at each particle's σ coordinate for fields on
half vertical levels (lower cell face, k+½), e.g. vertical velocity w. Layer k in the
field stores the value at the lower face of layer k (σ at k+½). The top face (σ=0) is
not stored and is implicitly zero by the kinematic boundary condition. Pass
`σ_levels_half` from `model.geometry`."""
function interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels, ::Face)
    (; npoints_output, north_pole_average, south_pole_average) = locator
    A_data = A.data
    arch = architecture(Aout)
    (; ring_starts, nlons, nlat) = geometry
    launch!(arch, LinearWorkOrder, (size(A_data, 2),), _compute_pole_averages_kernel!,
            north_pole_average, south_pole_average, A_data, ring_starts, nlons, nlat)
    launch!(
        arch, LinearWorkOrder, (npoints_output,),
        _interpolate_3D_face_kernel!,
        Aout, A_data, locator, positions, σ_levels,
    )
    return Aout
end
