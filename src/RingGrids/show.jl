# For grids, also add info about the number of rings, e.g.
# julia> A = randn(OctahedralClenshawGrid, 24)
# 3056-element, 47-ring OctahedralClenshawGrid{Float64}:
function Base.array_summary(io::IO, grid::AbstractGridArray, inds::Tuple{Vararg{Base.OneTo}})
    print(io, Base.dims2string(length.(inds)), ", $(get_nlat(grid))-ring ")
    Base.showarg(io, grid, true)
end