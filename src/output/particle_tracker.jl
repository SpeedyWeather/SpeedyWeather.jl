export ParticleTracker

"""
A ParticleTracker is implemented as a callback to output the trajectories
of particles from particle advection. To be added like

    add!(model.callbacks, ParticleTracker(spectral_grid; kwargs...))

Output done via netCDF. Fields and options are
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct ParticleTracker{NF} <: AbstractCallback
    "[OPTION] when to schedule particle tracking"
    schedule::Schedule = Schedule(every=Hour(3))

    "[OPTION] File name for netCDF file"
    file_name::String = "particles.nc"
    
    # COMPRESSION OPTIONS
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
    lon::Vector{NF} = zeros(NF, nparticles)
    lat::Vector{NF} = zeros(NF, nparticles)
    σ::Vector{NF} = zeros(NF, nparticles)
end

ParticleTracker(SG::SpectralGrid; kwargs...) =
    ParticleTracker{SG.NF}(;nparticles = SG.nparticles, kwargs...)

function initialize!(
    callback::ParticleTracker{NF},
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
) where NF

    initialize!(callback.schedule, progn.clock)

    # if model.output doesn't output create a folder anyway to store the particles.nc file
    if model.output.output == false
        (;output, feedback) = model
        output.id = get_run_id(output.path, output.id)
        output.run_path = create_output_folder(output.path, output.id) 
        
        feedback.id = output.id         # synchronize with feedback struct
        feedback.run_path = output.run_path
        feedback.progress_meter.desc = "Weather is speedy: run $(output.id) "
    end

    # CREATE NETCDF FILE, vector of NcVars for output
    (; run_path) = model.output
    (; file_name) = callback
    # model.output.output || @info "ParticleTracker stores data in $(joinpath(pwd(),file_name))"
    dataset = NCDataset(joinpath(run_path, file_name), "c")
    callback.netcdf_file = dataset
        
    # DEFINE NETCDF DIMENSIONS TIME
    (; time) = progn.clock
    time_string = "hours since $(Dates.format(time, "yyyy-mm-dd HH:MM:0.0"))"
    defDim(dataset, "time", Inf)       # unlimited time dimension
    defVar(dataset, "time", Float64, ("time",),
        attrib=Dict("units"=>time_string, "long_name"=>"time",
        "standard_name"=>"time", "calendar"=>"proleptic_gregorian"))

    # PARTICLE DIMENSION
    nparticles = length(progn.particles)
    callback.nparticles = nparticles
    defDim(dataset, "particle", nparticles)
    defVar(dataset, "particle", Int64, ("particle",),
        attrib=Dict("units"=>"1", "long_name"=>"particle identification number"))

    # coordinates of particles (the variables inside netCDF)
    defVar(dataset, "lon", NF, ("particle", "time"), attrib = 
        Dict("long_name"=>"longitude", "units"=>"degrees_north"),
        deflatelevel=callback.compression_level, shuffle=callback.shuffle)

    defVar(dataset, "lat", NF, ("particle", "time"), attrib = 
        Dict("long_name"=>"latitude", "units"=>"degrees_east"),
        deflatelevel=callback.compression_level, shuffle=callback.shuffle)

    defVar(dataset, "sigma", NF, ("particle", "time"), attrib = 
        Dict("long_name"=>"vertical sigma coordinate", "units"=>"1"),
        deflatelevel=callback.compression_level, shuffle=callback.shuffle)

    # pull particle locations into output work arrays
    for (p,particle) in enumerate(progn.particles)
        callback.lon[p] = particle.lon
        callback.lat[p] = particle.lat
        callback.σ[p]   = particle.σ
    end
    
    # rounding
    (; keepbits) = callback
    round!(callback.lon, keepbits)
    round!(callback.lat, keepbits)
    round!(callback.σ,   keepbits)
        
    # set initial conditions
    callback.netcdf_file["time"][1]     = 0.0
    callback.netcdf_file["lon"][:, 1]   = callback.lon
    callback.netcdf_file["lat"][:, 1]   = callback.lat
    callback.netcdf_file["sigma"][:, 1] = callback.σ
    NCDatasets.sync(callback.netcdf_file)
end

function callback!(
    callback::ParticleTracker,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    isscheduled(callback.schedule, progn.clock) || return nothing   # else escape immediately
    i = callback.schedule.counter+1     # +1 for initial conditions (not scheduled)

    # pull particle locations into output work arrays
    for (p,particle) in enumerate(progn.particles)
        callback.lon[p] = particle.lon
        callback.lat[p] = particle.lat
        callback.σ[p]   = particle.σ
    end
    
    # rounding
    (; keepbits) = callback
    round!(callback.lon, keepbits)
    round!(callback.lat, keepbits)
    round!(callback.σ,   keepbits)
        
    # write current particle locations to file
    (;time, start) = progn.clock
    time_passed_hrs = Millisecond(time - start).value/3600_000     # [ms] -> [hrs]
    callback.netcdf_file["time"][i]     = time_passed_hrs
    callback.netcdf_file["lon"][:, i]   = callback.lon
    callback.netcdf_file["lat"][:, i]   = callback.lat
    callback.netcdf_file["sigma"][:, i] = callback.σ
    NCDatasets.sync(callback.netcdf_file)
end

finish!(callback::ParticleTracker,args...) = NCDatasets.close(callback.netcdf_file)