export interpolate_3D!

# σ outside [σ_levels_full[1], σ_levels_full[end]] is pinned, not extrapolated
@inline function find_vertical_bracket(σ, σ_levels_full)
    n = length(σ_levels_full)
    k_lo = 1
    for k in 1:(n - 1)
        σ_levels_full[k] <= σ && (k_lo = k)
    end
    k_hi = k_lo + 1
    Δσ = σ_levels_full[k_hi] - σ_levels_full[k_lo]
    α = clamp((σ - σ_levels_full[k_lo]) / Δσ, zero(σ), one(σ))
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

@kernel inbounds = true function _interpolate_3D_kernel!(
        Aout, A_data,
        ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys,
        positions, σ_levels_full, north_pole_average, south_pole_average,
    )
    i = @index(Global, Linear)
    k_lo, k_hi, α = find_vertical_bracket(positions[i].σ, σ_levels_full)

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

"""$(TYPEDSIGNATURES)
Vertically blended anvil interpolation at each particle's σ coordinate.
Pole ring averages are stored in `locator.north_pole_average` and
`locator.south_pole_average` (pre-allocated via the 3D `AnvilLocator` constructor)
and reused across calls without CPU allocations."""
function interpolate_3D!(Aout, A, locator, geometry, positions, σ_levels_full)
    (; ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys, npoints_output) = locator
    (; north_pole_average, south_pole_average) = locator
    A_data = A.data
    arch = architecture(Aout)
    (; ring_starts, nlons, nlat) = geometry
    launch!(arch, LinearWorkOrder, (size(A_data, 2),), _compute_pole_averages_kernel!,
            north_pole_average, south_pole_average, A_data, ring_starts, nlons, nlat)
    launch!(
        arch, LinearWorkOrder, (npoints_output,),
        _interpolate_3D_kernel!,
        Aout, A_data, ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys,
        positions, σ_levels_full, north_pole_average, south_pole_average,
    )
    return Aout
end
