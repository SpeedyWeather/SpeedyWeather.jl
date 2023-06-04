# For grids, also add info about the number of rings, e.g.
# julia> A = randn(OctahedralClenshawGrid,24)
# 3056-element, 47-ring OctahedralClenshawGrid{Float64}:
function Base.array_summary(io::IO, grid::AbstractGrid, inds::Tuple{Vararg{Base.OneTo}})
    print(io, Base.dims2string(length.(inds)), ", $(get_nlat(grid))-ring ")
    Base.showarg(io, grid, true)
end

function plot(A::AbstractGrid;title::String="$(get_nlat(A))-ring $(typeof(A))")
    A_full = interpolate(full_grid(typeof(A)),A.nlat_half,A)
    plot(A_full;title)
end

function plot(A::AbstractFullGrid;title::String="$(get_nlat(A))-ring $(typeof(A))")

    A_matrix = Matrix(A)
    nlon,nlat = size(A_matrix)
    A_view = view(A_matrix,:,nlat:-1:1)

    plot_kwargs = pairs((   xlabel="˚E",
                            xfact=360/(nlon-1),
                            ylabel="˚N",
                            yfact=180/(nlat-1),
                            yoffset=-90,
                            title=title,
                            colormap=:viridis,
                            compact=true,
                            colorbar=true,
                            width=60,
                            height=30))

    UnicodePlots.heatmap(A_view';plot_kwargs...)
end