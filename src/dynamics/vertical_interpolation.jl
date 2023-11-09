"""
$(TYPEDSIGNATURES)
Interpolation weights for full to half level interpolation
on sigma coordinates. Following Fortran SPEEDY documentation eq. (1)."""
function σ_interpolation_weights(
    σ_levels_full::AbstractVector,
    σ_levels_half::AbstractVector)

    weights = zero(σ_levels_full)
    nlev = length(weights)

    for k in 1:nlev-1
        weights[k] = (log(σ_levels_half[k+1]) - log(σ_levels_full[k])) /
                        (log(σ_levels_full[k+1]) - log(σ_levels_full[k]))
    end
    # was log(0.99) in Fortran SPEEDY but doesn't make sense to me
    weights[end] =  (log(σ_levels_half[nlev+1]) - log(σ_levels_full[nlev])) /
                        (log(σ_levels_full[nlev]) - log(σ_levels_full[nlev-1]))
    
    return weights
end


"""
$(TYPEDSIGNATURES)
Given some generic column variable A defined at full levels, do a linear interpolation in
log(σ) to calculate its values at half-levels.
"""
function vertical_interpolate!(
    A_half::Vector,             # quantity A on half levels (excl top)
    A_full::Vector,             # quantity A on full levels
    G::Geometry,
)
    nlev = length(A_half)
    weights = G.full_to_half_interpolation

    # full levels contain one more for surface
    # TODO this is currently confusing because the surface fluxes use full[end]
    # as surface value which is technically on half levels though!
    @boundscheck nlev <= length(A_full) || throw(BoundsError)
    @boundscheck nlev <= length(weights) || throw(BoundsError)

    # For A at each full level k, compute A at the half-level below, i.e. at the boundary
    # between the full levels k and k+1. Fortran SPEEDY documentation eq. (1)
    for k = 1:nlev-1
        A_half[k] = A_full[k] + weights[k]*(A_full[k+1] - A_full[k])
    end

    # Compute the values at the surface separately
    A_half[nlev] = A_full[nlev] + weights[nlev]*(A_full[nlev] - A_full[nlev-1])

    return nothing
end