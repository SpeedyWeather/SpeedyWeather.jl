function plot(L::LowerTriangularMatrix{T}; mode::Function=abs) where T

    l,m = size(L)
    title ="$l×$m LowerTriangularMatrix{$T}"

    Lplot = similar(L,real(T))
    for lm in eachharmonic(L)
        Lplot[lm] = mode(L[lm]) 
    end

    # use at most 33x32 points in height x width, but fewer for smaller matrices
    height = min(l,33)
    width = min(m,32)

    plot_kwargs = pairs((   xlabel="m",
                            xoffset=-1,
                            ylabel="ℓ",
                            yoffset=-1,
                            title=title,
                            colormap=:inferno,
                            compact=true,
                            colorbar=true,
                            zlabel=string(mode),
                            array=true,
                            width=width,
                            height=height))

    return UnicodePlots.heatmap(Lplot;plot_kwargs...)
end