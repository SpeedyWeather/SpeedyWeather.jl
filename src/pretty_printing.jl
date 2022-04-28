function Base.show(io::IO, P::PrognosticVariables)

    lf = 1                      # first leapfrog index
    lev = size(P.temp)[end]     # surface level
    temp = view(P.temp,:,:,lf,lev)
    
    temp_surf_grid_degC = gridded(temp)     # to grid space
    temp_surf_grid_degC .-= 273.15          # and ˚C
    temp_surf_grid_degC = temp_surf_grid_degC[:,end:-1:1]       # flip latitudes

    nlon,nlat = size(temp_surf_grid_degC)

    plot_kwargs = pairs((   xlabel="˚E",
                            xfact=360/(nlon-1),
                            ylabel="˚N",
                            yfact=180/(nlat-1),
                            yoffset=-90,
                            zlabel="˚C",
                            title="Surface temperature",
                            colormap=:plasma))

    print(io,UnicodePlots.heatmap(temp_surf_grid_degC';plot_kwargs...))
end
    
    