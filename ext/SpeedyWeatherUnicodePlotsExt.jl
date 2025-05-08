module SpeedyWeatherUnicodePlotsExt

using SpeedyWeather, UnicodePlots
using DocStringExtensions

"""$(TYPEDSIGNATURES)
A UnicodePlots heatmap visualising the elements in a LowerTriangularMatrix.
Takes by default the `mode=abs` absolute value of every element to show the magnitude
of complex spherical harmonic coefficients."""
function UnicodePlots.heatmap(L::LowerTriangularMatrix{T}; mode::Function=abs) where T

    l, m = size(L, as=Matrix)
    title ="$l×$m LowerTriangularMatrix{$T}"

    Lplot = similar(L, real(T))
    for lm in eachharmonic(L)
        Lplot[lm] = mode(L[lm]) 
    end
    Lplot = Matrix(Lplot)

    # use at most 33x32 points in height x width, but fewer for smaller matrices
    height = min(l, 33)
    width = min(m, 32)

    plot_kwargs = pairs((   xlabel="degree m",
                            xoffset=-1,
                            ylabel="order l",
                            yoffset=-1,
                            title=title,
                            colormap=:inferno,
                            compact=true,
                            colorbar=true,
                            zlabel=string(mode),
                            array=true,
                            width=width,
                            height=height))

    return UnicodePlots.heatmap(Lplot; plot_kwargs...)
end

"""$(TYPEDSIGNATURES)
A UnicodePlots heatmap visualising the data on an `AbstractGrid`.
General method that interpolates (from a reduced grid) onto a full grid
so that it can be visualised as a matrix."""
function UnicodePlots.heatmap(A::RingGrids.AbstractGrid; title::String="$(get_nlat(A))-ring $(typeof(A))")
    A_full = RingGrids.interpolate(full_grid_type(A), A.nlat_half, A)
    plot(A_full; title)
end

"""$(TYPEDSIGNATURES)
A UnicodePlots heatmap visualising the data on an `AbstractFullGrid`."""
function UnicodePlots.heatmap(A::RingGrids.AbstractFullGrid; title::String="$(get_nlat(A))-ring $(typeof(A))")

    A_matrix = Matrix(A)
    nlon, nlat = size(A_matrix)
    A_view = view(A_matrix, :, nlat:-1:1)

    # use at most 30 points in height, but fewer for smaller grids
    # small grids are then displayed as 1 character per grid point
    height = min(nlat, 30)
    width = 2height

    plot_kwargs = pairs((   xlabel="˚E",
                            xfact=360/(nlon-1),
                            ylabel="˚N",
                            yfact=180/(nlat-1),
                            yoffset=-90,
                            title=title,
                            colormap=:viridis,
                            compact=true,
                            colorbar=true,
                            width=width,
                            height=height))

    return UnicodePlots.heatmap(A_view'; plot_kwargs...)
end
end # module