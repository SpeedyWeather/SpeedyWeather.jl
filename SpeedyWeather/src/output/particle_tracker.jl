export ParticleTracker

"""
A ParticleTracker is implemented as a callback to output the trajectories
of particles from particle advection. To be added like

    add!(model.callbacks, ParticleTracker(spectral_grid; kwargs...))

Output done via netCDF. Fields and options are
$(TYPEDFIELDS)"""
@kwdef mutable struct ParticleTracker{V} <: AbstractCallback
    "[OPTION] when to schedule particle tracking"
    schedule::Schedule = Schedule(every = Hour(4))

    "[OPTION] File name for netCDF file"
    filename::String = "particles.nc"

    "[OPTION] Path for netCDF file, uses model.output.run_path if not specified"
    path::String = ""

    "[OPTION] lossless compression level; 1=low but fast, 9=high but slow"
    compression_level::Int = 1

    "[OPTION] shuffle/bittranspose filter for compression"
    shuffle::Bool = false

    "[OPTION] mantissa bits to keep, (14, 15, 16) means at least (2km, 1km, 500m) accurate locations"
    keepbits::Int = 15

    "Number of particles to track"
    nparticles::Int = 0

    "The netcdf file to be written into, will be created at initialization"
    netcdf_file::Union{NCDataset, Nothing} = nothing

    # tracking arrays
    lon::V
    lat::V
    σ::V
end

function ParticleTracker(SG::SpectralGrid; kwargs...)
    # dummy vector
    v = zeros(SG.NF, 0)
    return ParticleTracker{SG.VectorType}(; lon = v, lat = v, σ = v, kwargs...)
end

function initialize!(
        particle_tracker::ParticleTracker,
        vars::Variables,
        model::AbstractModel,
    )

    (; clock) = vars.prognostic
    initialize!(particle_tracker.schedule, clock)

    # CREATE NETCDF FILE, vector of NcVars for output
    (; filename) = particle_tracker
    path = particle_tracker.path == "" ? model.output.run_path : particle_tracker.path
    mkpath(path)
    model.output.active || @info "ParticleTracker writes to $(joinpath(path, filename)) as output=false"

    dataset = NCDataset(joinpath(path, filename), "c")
    particle_tracker.netcdf_file = dataset

    # DEFINE NETCDF DIMENSIONS TIME
    (; time) = clock
    time_string = "hours since $(Dates.format(time, "yyyy-mm-dd HH:MM:0.0"))"
    defDim(dataset, "time", Inf)       # unlimited time dimension
    defVar(
        dataset, "time", Float64, ("time",),
        attrib = Dict(
            "units" => time_string, "long_name" => "time",
            "standard_name" => "time", "calendar" => "proleptic_gregorian"
        )
    )

    # allocate tracking arrays now of correct size
    arch = model.spectral_grid.architecture
    NF = model.spectral_grid.NF
    (; particles) = vars.prognostic
    nparticles = length(particles)
    particle_tracker.nparticles = nparticles
    particle_tracker.lon = on_architecture(arch, zeros(NF, nparticles))
    particle_tracker.lat = on_architecture(arch, zeros(NF, nparticles))
    particle_tracker.σ = on_architecture(arch, zeros(NF, nparticles))

    # PARTICLE DIMENSION
    defDim(dataset, "particle", nparticles)
    defVar(
        dataset, "particle", Int64, ("particle",),
        attrib = Dict("units" => "1", "long_name" => "particle identification number")
    )

    # coordinates of particles (the variables inside netCDF)
    defVar(
        dataset, "lon", NF, ("particle", "time"), attrib =
            Dict("long_name" => "longitude", "units" => "degrees_north"),
        deflatelevel = particle_tracker.compression_level, shuffle = particle_tracker.shuffle
    )

    defVar(
        dataset, "lat", NF, ("particle", "time"), attrib =
            Dict("long_name" => "latitude", "units" => "degrees_east"),
        deflatelevel = particle_tracker.compression_level, shuffle = particle_tracker.shuffle
    )

    defVar(
        dataset, "sigma", NF, ("particle", "time"), attrib =
            Dict("long_name" => "vertical sigma coordinate", "units" => "1"),
        deflatelevel = particle_tracker.compression_level, shuffle = particle_tracker.shuffle
    )

    # pull particle locations into output work arrays
    for (p, particle) in enumerate(particles)
        particle_tracker.lon[p] = particle.lon
        particle_tracker.lat[p] = particle.lat
        particle_tracker.σ[p] = particle.σ
    end

    # rounding
    (; keepbits) = particle_tracker
    round!(particle_tracker.lon, keepbits)
    round!(particle_tracker.lat, keepbits)
    round!(particle_tracker.σ, keepbits)

    # set initial conditions
    particle_tracker.netcdf_file["time"][1] = 0.0
    particle_tracker.netcdf_file["lon"][:, 1] = particle_tracker.lon
    particle_tracker.netcdf_file["lat"][:, 1] = particle_tracker.lat
    particle_tracker.netcdf_file["sigma"][:, 1] = particle_tracker.σ
    return NCDatasets.sync(particle_tracker.netcdf_file)
end

function callback!(
        particle_tracker::ParticleTracker,
        vars::Variables,
        model::AbstractModel,
    )
    isscheduled(particle_tracker.schedule, vars.prognostic.clock) || return nothing   # else escape immediately
    i = particle_tracker.schedule.counter + 1     # +1 for initial conditions (not scheduled)

    # pull particle locations into output work arrays
    for (p, particle) in enumerate(vars.prognostic.particles)
        particle_tracker.lon[p] = particle.lon
        particle_tracker.lat[p] = particle.lat
        particle_tracker.σ[p] = particle.σ
    end

    # rounding
    (; keepbits) = particle_tracker
    round!(particle_tracker.lon, keepbits)
    round!(particle_tracker.lat, keepbits)
    round!(particle_tracker.σ, keepbits)

    # write current particle locations to file
    (; time, start) = vars.prognostic.clock
    time_passed_hrs = Millisecond(time - start).value / 3600_000     # [ms] -> [hrs]
    particle_tracker.netcdf_file["time"][i] = time_passed_hrs
    particle_tracker.netcdf_file["lon"][:, i] = particle_tracker.lon
    particle_tracker.netcdf_file["lat"][:, i] = particle_tracker.lat
    particle_tracker.netcdf_file["sigma"][:, i] = particle_tracker.σ
    return NCDatasets.sync(particle_tracker.netcdf_file)
end

finalize!(particle_tracker::ParticleTracker, args...) = NCDatasets.close(particle_tracker.netcdf_file)
