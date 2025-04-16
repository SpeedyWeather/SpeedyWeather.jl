"""$(TYPEDSIGNATURES)
Vertices are defined for every grid point on a ring grid through 4 points: east, south, west, north.
    
    - east: longitude mid-point with the next grid point east
    - south: longitude mid-point between the two closest grid points on one ring to the south
    - west: longitude mid-point with the next grid point west
    - north: longitude mid-point between the two closest grid points on one ring to the north

Example

       o ----- n ------ o

    o --- w --- c --- e --- o

         o ----- s ------ o

with cell center c (the grid point), e, s, w, n the vertices and o the surrounding grid points.
Returns 2xnpoints arrays for east, south, west, north each containing the longitude and latitude of the vertices."""
function get_vertices(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)

    npoints = get_npoints2D(Grid, nlat_half)
    londs, latds = get_londlatds(Grid, nlat_half)

    # use the locator to find the neighbouring grid points
    # to control the direction of "neighbouring" use a small offset Δ
    # to add/subtract from the coordinates
    Δ = 1e-4
    
    I = AnvilInterpolator(Grid, nlat_half, npoints)
    update_locator!(I, londs .+ Δ, latds .+ Δ)

    east = zeros(2, npoints)
    south = zeros(2, npoints)
    west = zeros(2, npoints)
    north = zeros(2, npoints)

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

    update_locator!(I, londs .- Δ, latds .- Δ)

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

    return east, south, west, north
end

"""$(TYPEDSIGNATURES)
Vertices for full grids, other definition than for reduced grids to prevent
a diamond shape of the cells. Use default rectangular instead.
Effectively rotating the vertices clockwise by 45˚, making east south-east etc."""
function get_vertices(Grid::Type{<:AbstractFullGridArray}, nlat_half::Integer)

    npoints = get_npoints2D(Grid, nlat_half)
    nlat = get_nlat(Grid, nlat_half)
    nlon = get_nlon(Grid, nlat_half)
    latd = get_latd(Grid, nlat_half)
    lond = get_lond(Grid, nlat_half)
    dλ = lond[2] - lond[1]

    # vertices for full grids are at north west, north east etc
    seast = zeros(2, npoints)
    swest = zeros(2, npoints)
    nwest = zeros(2, npoints)
    neast = zeros(2, npoints)

    # longitude is just shifted
    west = mod.(lond .- dλ/2, 360)
    east = mod.(lond .+ dλ/2, 360)

    @inbounds for j in 1:nlat
        ring = nlon*(j-1) + 1 : nlon*j
        nwest[1, ring] = west
        swest[1, ring] = west
        neast[1, ring] = east
        seast[1, ring] = east

        # average ring latitudes
        φ_north = j == 1 ? 90 : (latd[j] + latd[j-1])/2
        φ_south = j == nlat ? -90 : (latd[j] + latd[j+1])/2
        nwest[2, ring] .= φ_north
        swest[2, ring] .= φ_south
        neast[2, ring] .= φ_north
        seast[2, ring] .= φ_south
    end

    # retain clockwise order
    return seast, swest, nwest, neast
end

"""$(TYPEDSIGNATURES)
Return a 5xN matrix `polygons` (or grid cells) of `NTuple{2, Float64}` where the first 4 rows
are the vertices (E, S, W, N) of every grid points ij in 1:N, row 5 is duplicated north vertex
to close the grid cell. Use keyword arguemnt `add_nan=true` (default `false`) to add a 6th row
with (NaN, NaN) to separate grid cells when drawing them as a continuous line with `vec(polygons)`."""
function get_gridcell_polygons(
    Grid::Type{<:AbstractGridArray},
    nlat_half::Integer;
    add_nan::Bool = false,
)
    npoints = get_npoints2D(Grid, nlat_half)

    # vertex east, south, west, north (i.e. clockwise for every grid point)
    E, S, W, N = get_vertices(Grid, nlat_half)

    @boundscheck size(N) == size(W) == size(S) == size(E) || throw(BoundsError("Vertices must have the same size"))
    @boundscheck size(N) == (2, npoints) || throw(BoundsError("Number of vertices and npoints do not agree"))

    # number of vertices = 4, 5 to close the polygon, 6 to add a nan
    # to prevent grid lines to be drawn between cells
    nvertices = add_nan ? 6 : 5

    # allocate polygons as tuple so that no data copy has to be made in Makie
    polygons = Matrix{NTuple{2, Float64}}(undef, nvertices, npoints)

    @inbounds for ij in 1:npoints
        polygons[1, ij] = (E[1, ij], E[2, ij])  # clockwise
        polygons[2, ij] = (S[1, ij], S[2, ij])
        polygons[3, ij] = (W[1, ij], W[2, ij])
        polygons[4, ij] = (N[1, ij], N[2, ij])
        polygons[5, ij] = (E[1, ij], E[2, ij])  # back to east to close the polygon        
    end

    if add_nan  # add a NaN to separate grid cells
        for ij in 1:npoints
            polygons[6, ij] = (NaN, NaN)
        end
    end

    return polygons
end

get_gridcell_polygons(grid::AbstractGridArray; kwargs...) =
    get_gridcell_polygons(typeof(grid), grid.nlat_half; kwargs...)