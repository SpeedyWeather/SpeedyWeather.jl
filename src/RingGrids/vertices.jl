"""$(TYPEDSIGNATURES)
Vertices are defined for every grid point on a ring grid through 4 points: north, west, south, east.
    
    - north: longitude mid-point between the two closest grid points on one ring to the north
    - east: longitude mid-point with the next grid point east
    - south: longitude mid-point between the two closest grid points on one ring to the south
    - west: longitude mid-point with the next grid point west

Example

       o ----- n ------ o

    o --- w --- c --- e --- o

         o ----- s ------ o

with cell center c (the grid point), n, e, s, w the vertices and o the surrounding grid points.
Returns 2xnpoints arrays for north, west, south, east, each containing the longitude and latitude of the vertices."""
function get_vertices(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)

    npoints = get_npoints2D(Grid, nlat_half)
    latds, londs = get_latdlonds(Grid, nlat_half)

    # use the locator to find the neighbouring grid points
    # to control the direction of "neighbouring" use a small offset Δ
    # to add/subtract from the coordinates
    Δ = 1e-4
    
    I = AnvilInterpolator(Grid, nlat_half, npoints)
    update_locator!(I, latds .+ Δ, londs .+ Δ)

    north = zeros(2, npoints)
    east = zeros(2, npoints)
    south = zeros(2, npoints)
    west = zeros(2, npoints)

    (; ij_as, ij_bs, ij_cs, ij_ds) = I.locator
    @inbounds for (ij, (a, b, c, d)) in enumerate(zip(ij_as, ij_bs, ij_cs, ij_ds))
        if a == 0                       # north pole but
            north[1, ij] = londs[ij]    # pick same longitude
            north[2, ij] = 90
        else
            # account for averages across prime meridian in 0..360˚ coordinates
            λa = londs[b] < londs[a] ? londs[a] - 360 : londs[a]
            λb = londs[b]
            north[1, ij] = mod((λa + λb)/2, 360)
            north[2, ij] = (latds[a] + latds[b])/2
        end

        # account for averages across prime meridian in 0..360˚ coordinates
        λc = londs[d] < londs[c] ? londs[c] - 360 : londs[c]
        λd = londs[d]
        east[1, ij] = mod((λc + λd)/2, 360)
        east[2, ij] = (latds[c] + latds[d])/2
    end

    update_locator!(I, latds .- Δ, londs .- Δ)

    (; ij_as, ij_bs, ij_cs, ij_ds) = I.locator
    @inbounds for (ij, (a, b, c, d)) in enumerate(zip(ij_as, ij_bs, ij_cs, ij_ds))
        if c == -1                    # south pole
            south[1, ij] = londs[ij]  # pick same longitude
            south[2, ij] = -90
        else
            # account for averages across prime meridian in 0..360˚ coordinates
            λc = londs[c] > londs[d] ? londs[c] - 360 : londs[c]
            λd = londs[d]
            south[1, ij] = mod((λc + λd)/2, 360)
            south[2, ij] = (latds[c] + latds[d])/2
        end

        # account for averages across prime meridian in 0..360˚ coordinates
        λa = londs[b] < londs[a] ? londs[a] - 360 : londs[a]
        λb = londs[b]
        west[1, ij] = mod((λa + λb)/2, 360)
        west[2, ij] = (latds[a] + latds[b])/2
    end

    return north, west, south, east
end