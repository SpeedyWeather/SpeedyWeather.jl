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

# adapted from ProgressMeter.jl
function speedstring(sec_per_iter,dt_in_sec)
    if sec_per_iter == Inf
        return "  N/A  days/day"
    end

    sim_time_per_time = dt_in_sec/sec_per_iter

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

# adapted from ProgressMeter.jl
function remaining_time(p::ProgressMeter.Progress)
    elapsed_time = time() - p.tinit
    est_total_time = elapsed_time * (p.n - p.start) / (p.counter - p.start)
    if 0 <= est_total_time <= typemax(Int)
        eta_sec = round(Int, est_total_time - elapsed_time )
        eta = ProgressMeter.durationstring(eta_sec)
    else
        eta = "N/A"
    end
    return eta
end

# hack: define global constant whose element will be changed in initialize_feedback
# used to pass on the time step to ProgressMeter.speedstring via calling this
# constant from the ProgressMeter module
const DT_IN_SEC = Ref(1800)

function ProgressMeter.speedstring(sec_per_iter,dt_in_sec=SpeedyWeather.DT_IN_SEC)
    speedstring(sec_per_iter,dt_in_sec[])
end

"""
    readable_secs(secs::Real) -> Dates.CompoundPeriod

Returns `Dates.CompoundPeriod` rounding to either (days, hours), (hours, minutes), (minutes,
seconds), or seconds with 1 decimal place accuracy for >10s and two for less.
E.g.
```julia
julia> readable_secs(12345)
3 hours, 26 minutes
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(secs * 10 ^ 3))
    if millisecs >= Dates.Day(1)
        return Dates.canonicalize(round(millisecs, Dates.Hour))
    elseif millisecs >= Dates.Hour(1)
        return Dates.canonicalize(round(millisecs, Dates.Minute))
    elseif millisecs >= Dates.Minute(1)
        return Dates.canonicalize(round(millisecs, Dates.Second))
    elseif millisecs >= Dates.Second(10)
        return Dates.canonicalize(round(millisecs, Dates.Millisecond(100)))
    end
    return Dates.canonicalize(round(millisecs, Dates.Millisecond(10)))
end