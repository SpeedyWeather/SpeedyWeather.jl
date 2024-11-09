import GeoMakie: GeometryBasics

function GeometryBasics.earcut_triangulate(polygon::Vector{Vector{Point{3,Float64}}})
    # Here, we are assuming that the polygon is actually planar,
    # but the plane is in 3D, and not necessarily the XY plane.
    # So, we can find the plane of best fit using the first three points of the polygon!
    p1, p2, p3 = polygon[1][1], polygon[1][2], polygon[1][3]
    v1 = p2 - p1
    v2 = p3 - p1
    normal = cross(v1, v2)
    x = v1
    y = cross(normal, x)

    projected_polygon = map(ring -> map(p -> Point2{Float64}(dot(p, x), dot(p, y)), ring), polygon)
    lengths = map(x -> UInt32(length(x)), projected_polygon)
    len = UInt32(length(lengths))
    array = ccall((:u32_triangulate_f64, GeometryBasics.libearcut), Tuple{Ptr{GeometryBasics.GLTriangleFace},Cint},
                  (Ptr{Ptr{Float64}}, Ptr{UInt32}, UInt32), projected_polygon, lengths, len)
    return unsafe_wrap(Vector{GeometryBasics.GLTriangleFace}, array[1], array[2]; own=true)
end

"""$(TYPEDSIGNATURES)
Transpose (and copy) the 4 vertices of every grid point to obtain the faces of the grid.
Return a 6xN matrix `faces` of Point2{Float64} where the first 4 rows are the vertices (N, E, S, W)
of every grid points ij in 1:N, row 5 is duplicated north vertex to close the grid cell,
row 6 is NaN to separate grid cells when drawing them as a continuous line with `vec(faces)`."""
function get_faces(
    Grid::Type{<:AbstractGridArray},
    nlat_half::Integer;
    add_nan::Bool = false,
)
    npoints = RingGrids.get_npoints2D(Grid, nlat_half)

    # vertex north, west, south, east (i.e. clockwise for every grid point)
    N, W, S, E = RingGrids.get_vertices(Grid, nlat_half)

    @boundscheck size(N) == size(W) == size(S) == size(E) || throw(BoundsError("Vertices must have the same size"))
    @boundscheck size(N) == (2, npoints) || throw(BoundsError("Number of vertices and npoints do not agree"))

    # number of vertices = 4, 5 to close the polygon, 6 to add a nan
    # to prevent grid lines to be drawn between cells
    nvertices = add_nan ? 6 : 5

    # allocate faces as Point2{Float64} so that no data copy has to be made in Makie
    faces = Matrix{Point2{Float64}}(undef, nvertices, npoints)

    @inbounds for ij in 1:npoints
        faces[1, ij] = Point2{Float64}(N[1, ij], N[2, ij])  # clockwise
        faces[2, ij] = Point2{Float64}(E[1, ij], E[2, ij])
        faces[3, ij] = Point2{Float64}(S[1, ij], S[2, ij])
        faces[4, ij] = Point2{Float64}(W[1, ij], W[2, ij])
        faces[5, ij] = Point2{Float64}(N[1, ij], N[2, ij])  # back to north to close the polygon        
    end

    if add_nan  # add a NaN to separate grid cells
        for ij in 1:npoints
            faces[6, ij] = Point2{Float64}(NaN, NaN)
        end
    end

    return faces
end

get_faces(grid::AbstractGridArray; kwargs...) = get_faces(typeof(grid), grid.nlat_half; kwargs...)