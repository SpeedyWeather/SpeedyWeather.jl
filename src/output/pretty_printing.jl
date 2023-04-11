function Base.show(io::IO, P::PrognosticVariables)

    ζ = P.layers[end].timesteps[1].vor   # create a view on vorticity
    ζ_grid = Matrix(gridded(ζ))         # to grid space
    ζ_grid = ζ_grid[:,end:-1:1]         # flip latitudes

    nlon,nlat = size(ζ_grid)

    plot_kwargs = pairs((   xlabel="˚E",
                            xfact=360/(nlon-1),
                            ylabel="˚N",
                            yfact=180/(nlat-1),
                            yoffset=-90,
                            title="Surface relative vorticity",
                            colormap=:viridis,
                            compact=true,
                            colorbar=true,
                            width=60,
                            height=30))

    print(io,UnicodePlots.heatmap(ζ_grid';plot_kwargs...))
end

# hack: define global constant whose element will be changed in initialize_feedback
# used to pass on the time step to ProgressMeter.speedstring via calling this
# constant from the ProgressMeter module
const DT_IN_SEC = [1800]

function ProgressMeter.speedstring(sec_per_iter,dt_in_sec=SpeedyWeather.DT_IN_SEC)
    if sec_per_iter == Inf
        return "  N/A  days/day"
    end

    sim_time_per_time = dt_in_sec[1]/sec_per_iter

    for (divideby, unit) in (   (365*1_000, "millenia"),
                                (365, "years"),
                                (1, "days"),
                                (1/24, "hours"))    
        if (sim_time_per_time / divideby) > 2
            return @sprintf "%5.2f %2s/day" (sim_time_per_time / divideby) unit
        end
    end
    return " <2 hours/days"
end