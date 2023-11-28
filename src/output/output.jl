"""
Number of mantissa bits to keep for each prognostic variable when compressed for
netCDF and .jld2 data output.
$(TYPEDFIELDS)"""
Base.@kwdef struct Keepbits
    u::Int = 7
    v::Int = 7
    vor::Int = 5
    div::Int = 5
    temp::Int = 10
    pres::Int = 12
    humid::Int = 7
    precip_cond::Int = 7
    precip_conv::Int = 7
    cloud::Int = 7
end

function Base.show(io::IO,K::Keepbits)
    println(io,"$(typeof(K))")
    keys = propertynames(K)
    print_fields(io,K,keys)
end

# default number format for output
const DEFAULT_OUTPUT_NF = Float32
const DEFAULT_OUTPUT_DT = Hour(6)

"""
$(TYPEDSIGNATURES)
NetCDF output writer. Contains all output options and auxiliary fields for output interpolation.
To be initialised with `OutputWriter(::SpectralGrid,::Type{<:ModelSetup},kwargs...)` to pass on the
resolution information and the model type which chooses which variables to output. Options include
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct OutputWriter{NF<:Union{Float32,Float64},Model<:ModelSetup} <: AbstractOutputWriter

    spectral_grid::SpectralGrid

    # FILE OPTIONS
    output::Bool = false                    
    
    "[OPTION] path to output folder, run_???? will be created within"
    path::String = pwd()
    
    "[OPTION] run identification number/string"
    id::String = "0001"
    run_path::String = ""                   # will be determined in initalize!
    
    "[OPTION] name of the output netcdf file"
    filename::String = "output.nc"
    
    "[OPTION] also write restart file if output==true?"
    write_restart::Bool = true
    pkg_version::VersionNumber = pkgversion(SpeedyWeather)

    # WHAT/WHEN OPTIONS
    startdate::DateTime = DateTime(2000,1,1)

    "[OPTION] output frequency, time step"
    output_dt::Second = DEFAULT_OUTPUT_DT

    "[OPTION] which variables to output, u, v, vor, div, pres, temp, humid"
    output_vars::Vector{Symbol} = default_output_vars(Model)

    "[OPTION] missing value to be used in netcdf output"
    missing_value::NF = NaN

    # COMPRESSION OPTIONS
    "[OPTION] lossless compression level; 1=low but fast, 9=high but slow"
    compression_level::Int = 3

    "[OPTION] shuffle/bittranspose filter for compression"
    shuffle::Bool = true
    
    "[OPTION] mantissa bits to keep for every variable"
    keepbits::Keepbits = Keepbits()

    # TIME STEPS AND COUNTERS (initialize later)
    output_every_n_steps::Int = 0           # output frequency
    timestep_counter::Int = 0               # time step counter
    output_counter::Int = 0                 # output step counter
    
    # the netcdf file to be written into, will be create
    netcdf_file::Union{NCDataset,Nothing} = nothing

    # INPUT GRID (the one used in the dynamical core)
    input_Grid::Type{<:AbstractGrid} = spectral_grid.Grid
    
    # Output as matrix (particularly for reduced grids)
    "[OPTION] sort grid points into a matrix (interpolation-free), for OctahedralClenshawGrid, OctaHEALPixGrid only"
    as_matrix::Bool = false

    quadrant_rotation::NTuple{4,Int} = (0,1,2,3)    # rotation of output quadrant
                                                    # matrix of output quadrant
    matrix_quadrant::NTuple{4,Tuple{Int,Int}} = ((2,2),(1,2),(1,1),(2,1))
    
    # OUTPUT GRID
    "[OPTION] the grid used for output, full grids only"
    output_Grid::Type{<:AbstractFullGrid} = RingGrids.full_grid(input_Grid)
    
    "[OPTION] the resolution of the output grid, default: same nlat_half as in the dynamical core"
    nlat_half::Int = spectral_grid.nlat_half
    nlon::Int = as_matrix ? RingGrids.matrix_size(input_Grid,spectral_grid.nlat_half)[1] :
                                                    RingGrids.get_nlon(output_Grid,nlat_half)
    nlat::Int =  as_matrix ? RingGrids.matrix_size(input_Grid,spectral_grid.nlat_half)[2] :
                                                    RingGrids.get_nlat(output_Grid,nlat_half)
    npoints::Int = nlon*nlat
    nlev::Int = spectral_grid.nlev
    interpolator::AbstractInterpolator = DEFAULT_INTERPOLATOR(NF,input_Grid,spectral_grid.nlat_half,npoints)

    # fields to output (only one layer, reuse over layers)
    const u::Matrix{NF} = fill(missing_value,nlon,nlat)
    const v::Matrix{NF} = fill(missing_value,nlon,nlat)
    const vor::Matrix{NF} = fill(missing_value,nlon,nlat)
    const div::Matrix{NF} = fill(missing_value,nlon,nlat)
    const temp::Matrix{NF} = fill(missing_value,nlon,nlat)
    const pres::Matrix{NF} = fill(missing_value,nlon,nlat)
    const humid::Matrix{NF} = fill(missing_value,nlon,nlat)
    const precip_cond::Matrix{NF} = fill(missing_value,nlon,nlat)
    const precip_conv::Matrix{NF} = fill(missing_value,nlon,nlat)
    const cloud::Matrix{NF} = fill(missing_value,nlon,nlat)
end

# generator function pulling grid resolution and time stepping from ::SpectralGrid and ::TimeStepper
function OutputWriter(
    spectral_grid::SpectralGrid,
    ::Type{Model};
    NF::Type{<:Union{Float32,Float64}} = DEFAULT_OUTPUT_NF,
    kwargs...
) where {Model<:ModelSetup}
    return OutputWriter{NF,Model}(;spectral_grid,kwargs...)
end

# default variables to output by model
default_output_vars(::Type{<:Barotropic}) = [:vor,:u,:v]
default_output_vars(::Type{<:ShallowWater}) = [:vor,:u,:v]
default_output_vars(::Type{<:PrimitiveDry}) = [:vor,:u,:v,:temp,:pres]
default_output_vars(::Type{<:PrimitiveWet}) = [:vor,:u,:v,:temp,:humid,:pres,:precip,:cloud]

# print all fields with type <: Number
function Base.show(io::IO,O::AbstractOutputWriter)
    println(io,"$(typeof(O))")
    keys = propertynames(O)

    # remove interpolator from being printed, TODO: implement show for AbstractInterpolator
    keys_filtered = filter(key -> ~(getfield(O,key) isa AbstractInterpolator),keys)
    print_fields(io,O,keys_filtered)
end

"""
$(TYPEDSIGNATURES)
Creates a netcdf file on disk and the corresponding `netcdf_file` object preallocated with output variables
and dimensions. `write_output!` then writes consecuitive time steps into this file.
"""
function initialize!(   
    output::OutputWriter{output_NF,Model},
    feedback::AbstractFeedback,
    time_stepping::TimeStepper,
    clock::Clock,
    diagn::DiagnosticVariables,
    model::Model
) where {output_NF,Model}
    
    output.output || return nothing     # exit immediately for no output
    
    # GET RUN ID, CREATE FOLDER
    # get new id only if not already specified
    output.id = get_run_id(output.path,output.id)
    output.run_path = create_output_folder(output.path,output.id) 
    
    feedback.id = output.id         # synchronize with feedback struct
    feedback.run_path = output.run_path
    feedback.progress_meter.desc = "Weather is speedy: run $(output.id) "
    feedback.output = true          # if output=true set feedback.output=true too!

    # OUTPUT FREQUENCY
    output.output_every_n_steps = max(1,round(Int,
            Millisecond(output.output_dt).value/time_stepping.Δt_millisec.value))
    output.output_dt = Second(round(Int,output.output_every_n_steps*time_stepping.Δt_sec))

    # RESET COUNTERS
    output.timestep_counter = 0         # time step counter
    output.output_counter = 0           # output step counter

    # CREATE NETCDF FILE, vector of NcVars for output
    (; run_path, filename, output_vars) = output
    dataset = NCDataset(joinpath(run_path,filename),"c")
    output.netcdf_file = dataset
    
    # DEFINE NETCDF DIMENSIONS TIME
    (;startdate) = output
    time_string = "hours since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    defDim(dataset,"time",Inf)       # unlimited time dimension
    defVar(dataset,"time",Float64,("time",),attrib=Dict("units"=>time_string,"long_name"=>"time",
            "standard_name"=>"time","calendar"=>"proleptic_gregorian"))
    
    # DEFINE NETCDF DIMENSIONS SPACE
    (;input_Grid, output_Grid, nlat_half) = output
    
    if output.as_matrix == false        # interpolate onto (possibly different) output grid
        lond = get_lond(output_Grid,nlat_half)
        latd = get_latd(output_Grid,nlat_half)
        nlon = length(lond)
        nlat = length(latd)
        lon_name, lon_units, lon_longname = "lon","degrees_east","longitude"
        lat_name, lat_units, lat_longname = "lat","degrees_north","latitude"
        
    else                                # output grid directly into a matrix (resort grid points, no interpolation)
        (;nlat_half) = diagn            # don't use output.nlat_half as not supported for output_matrix
        nlon,nlat = RingGrids.matrix_size(input_Grid,nlat_half)     # size of the matrix output
        lond = collect(1:nlon)                                      # just enumerate grid points for lond, latd
        latd = collect(1:nlat)
        lon_name, lon_units, lon_longname = "i","1","horizontal index i"
        lat_name, lat_units, lat_longname = "j","1","horizontal index j"
    end
    
    # INTERPOLATION: PRECOMPUTE LOCATION INDICES
    latds, londs = RingGrids.get_latdlonds(output_Grid,output.nlat_half)
    output.as_matrix || RingGrids.update_locator!(output.interpolator,latds,londs)
        
    σ = output_NF.(model.geometry.σ_levels_full)
    defVar(dataset,lon_name,lond,(lon_name,),attrib=Dict("units"=>lon_units,"long_name"=>lon_longname))
    defVar(dataset,lat_name,latd,(lat_name,),attrib=Dict("units"=>lat_units,"long_name"=>lat_longname))
    defVar(dataset,"lev",σ,("lev",),attrib=Dict("units"=>"1","long_name"=>"sigma levels"))
    
    # VARIABLES, define every variable here that could be output
    (;compression_level) = output
    missing_value = convert(output_NF,output.missing_value)
    
    # given pres the right name, depending on ShallowWaterModel or PrimitiveEquationModel
    pres_name = Model <: ShallowWater ? "interface displacement" : "surface pressure"
    pres_unit = Model <: ShallowWater ? "m" : "hPa"
    
    # zonal wind
    u_attribs = Dict("long_name"=>"zonal wind","units"=>"m/s","_FillValue"=>missing_value)
    :u in output_vars && defVar(dataset,"u",output_NF,(lon_name,lat_name,"lev","time"),attrib=u_attribs,
                                deflatelevel=compression_level,shuffle=output.shuffle)

    # meridional wind
    v_attribs = Dict("long_name"=>"meridional wind","units"=>"m/s","_FillValue"=>missing_value)
    :v in output_vars && defVar(dataset,"v",output_NF,(lon_name,lat_name,"lev","time"),attrib=v_attribs,
                                deflatelevel=compression_level,shuffle=output.shuffle)

    # vorticity
    vor_attribs = Dict("long_name"=>"relative vorticity","units"=>"1/s","_FillValue"=>missing_value)
    :vor in output_vars && defVar(dataset,"vor",output_NF,(lon_name,lat_name,"lev","time"),attrib=vor_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)

    # divergence
    div_attribs = Dict("long_name"=>"divergence","units"=>"1/s","_FillValue"=>missing_value)
    :div in output_vars && defVar(dataset,"div",output_NF,(lon_name,lat_name,"lev","time"),attrib=div_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)

    # pressure / interface displacement
    pres_attribs = Dict("long_name"=>pres_name,"units"=>pres_unit,"_FillValue"=>missing_value)
    :pres in output_vars && defVar(dataset,"pres",output_NF,(lon_name,lat_name,"time"),attrib=pres_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)

    # temperature
    temp_attribs = Dict("long_name"=>"temperature","units"=>"degC","_FillValue"=>missing_value)
    :temp in output_vars && defVar(dataset,"temp",output_NF,(lon_name,lat_name,"lev","time"),attrib=temp_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)
                                    
    # humidity
    humid_attribs = Dict("long_name"=>"specific humidity","units"=>"kg/kg","_FillValue"=>missing_value)
    :humid in output_vars && defVar(dataset,"humid",output_NF,(lon_name,lat_name,"lev","time"),attrib=humid_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)

    # orography
    if :orography in output_vars    # write orography directly to file
        orog_attribs = Dict("long_name"=>"orography","units"=>"m","_FillValue"=>missing_value)
        orog_grid = model.orography.orography
        orog_matrix = output.u
        output.as_matrix && (orog_matrix = Matrix(orog_grid))
        output.as_matrix || RingGrids.interpolate!(output_Grid(output.u),orog_grid,output.interpolator)
        defVar(dataset,"orography",orog_matrix,(lon_name,lat_name),attrib=orog_attribs)
    end

    # large-scale condensation
    precip_cond_attribs = Dict("long_name"=>"large-scale precipitation","units"=>"mm/dt","_FillValue"=>missing_value)
    :precip in output_vars && defVar(dataset,"precip_cond",output_NF,(lon_name,lat_name,"time"),attrib=precip_cond_attribs,
                                        deflatelevel=compression_level,shuffle=output.shuffle)

    # convective precipitation
    precip_conv_attribs = Dict("long_name"=>"convective precipitation","units"=>"mm/dt","_FillValue"=>missing_value)
    :precip in output_vars && defVar(dataset,"precip_conv",output_NF,(lon_name,lat_name,"time"),attrib=precip_conv_attribs,
                                        deflatelevel=compression_level,shuffle=output.shuffle)
    
    # convective precipitation
    cloud_top_attribs = Dict("long_name"=>"cloud top","units"=>"Pa","_FillValue"=>missing_value)
    :cloud in output_vars && defVar(dataset,"cloud_top",output_NF,(lon_name,lat_name,"time"),attrib=cloud_top_attribs,
                                    deflatelevel=compression_level,shuffle=output.shuffle)

    # WRITE INITIAL CONDITIONS TO FILE
    write_netcdf_variables!(output,diagn)
    write_netcdf_time!(output,clock.time)

    # also export parameters into run????/parameters.txt
    parameters_txt = open(joinpath(output.run_path,"parameters.txt"),"w")
    println(parameters_txt,model.spectral_grid)
    println(parameters_txt,model.planet)
    println(parameters_txt,model.atmosphere)
    println(parameters_txt,model.time_stepping)
    println(parameters_txt,model.output)
    println(parameters_txt,model.initial_conditions)
    println(parameters_txt,model.horizontal_diffusion)
    model isa Union{ShallowWater,PrimitiveEquation} && println(parameters_txt,model.implicit)
    model isa Union{ShallowWater,PrimitiveEquation} && println(parameters_txt,model.orography)
    close(parameters_txt)
end

"""
$(TYPEDSIGNATURES)
Checks existing `run_????` folders in `path` to determine a 4-digit `id` number
by counting up. E.g. if folder run_0001 exists it will return the string "0002".
Does not create a folder for the returned run id.
"""
function get_run_id(path::String,id::String)
    # if run_???? folder doesn't exist yet don't change the id
    run_id = string("run_",run_id_to_string(id))
    !(run_id in readdir(path)) && return id

    # otherwise pull list of existing run_???? folders via readdir
    pattern = r"run_\d\d\d\d"               # run_???? in regex
    runlist = filter(x->startswith(x,pattern),readdir(path))
    runlist = filter(x->endswith(  x,pattern),runlist)
    existing_runs = [parse(Int,id[5:end]) for id in runlist]

    # get the run id from existing folders
    if length(existing_runs) == 0           # if no runfolder exists yet
        run_id = 1                          # start with run_0001
    else
        run_id = maximum(existing_runs)+1   # next run gets id +1
    end
    
    return @sprintf("%04d",run_id)
end

"""
$(TYPEDSIGNATURES)
Creates a new folder `run_*` with the identification `id`. Also returns the full path
`run_path` of that folder.
"""
function create_output_folder(path::String,id::Union{String,Int})
    run_id = string("run_",run_id_to_string(id))
    run_path = joinpath(path,run_id)
    @assert !(run_id in readdir(path)) "Run folder $run_path already exists."
    mkdir(run_path)             # actually create the folder
    return run_path
end

run_id_to_string(run_id::Integer) = @sprintf("%04d",run_id)
run_id_to_string(run_id::String) = run_id


"""
$(TYPEDSIGNATURES)
Writes the variables from `diagn` of time step `i` at time `time` into `outputter.netcdf_file`.
Simply escapes for no netcdf output of if output shouldn't be written on this time step.
Interpolates onto output grid and resolution as specified in `outputter`, converts to output
number format, truncates the mantissa for higher compression and applies lossless compression."""
function write_output!( outputter::OutputWriter,        # everything for netcdf output
                        time::DateTime,                 # model time for output
                        diagn::DiagnosticVariables)     # all diagnostic variables

    outputter.timestep_counter += 1                                 # increase counter
    (; output, output_every_n_steps, timestep_counter ) = outputter
    output || return nothing                                        # escape immediately for no netcdf output
    timestep_counter % output_every_n_steps == 0 || return nothing  # escape if output not written on this step

    # WRITE VARIABLES
    write_netcdf_variables!(outputter,diagn)
    write_netcdf_time!(outputter,time)
end

"""
$(TYPEDSIGNATURES)
Write the current time `time::DateTime` to the netCDF file in `output::OutputWriter`."""
function write_netcdf_time!(output::OutputWriter,
                            time::DateTime)
    
    (; netcdf_file, startdate ) = output
    i = output.output_counter

    time_passed = Millisecond(time-startdate)
    time_hrs = time_passed.value/3600_000       # [ms] -> [hrs]
    netcdf_file["time"][i] = time_hrs
    NCDatasets.sync(netcdf_file)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Write diagnostic variables from `diagn` to the netCDF file in `output::OutputWriter`."""
function write_netcdf_variables!(   output::OutputWriter,
                                    diagn::DiagnosticVariables{NF,Grid,Model}) where {NF,Grid,Model}

    output.output_counter += 1                  # increase output step counter
    (;output_vars) = output                     # Vector{Symbol} of variables to output
    i = output.output_counter

    (;u, v, vor, div, pres, temp, humid, precip_cond, precip_conv, cloud) = output
    (;output_Grid, interpolator) = output
    (;quadrant_rotation, matrix_quadrant) = output
    (;netcdf_file, keepbits) = output

    for (k,diagn_layer) in enumerate(diagn.layers)
        
        (; u_grid, v_grid, vor_grid, div_grid, temp_grid, humid_grid ) = diagn_layer.grid_variables

        if output.as_matrix     # resort gridded variables interpolation-free into a matrix

            # create (matrix,grid) tuples for simultaneous grid -> matrix conversion
            # TODO this currently does the Matrix! conversion to all variables, not just output_vars
            # as arrays are always initialised  
            MGs = ((M,G) for (M,G) in zip((u,v,vor,div,temp,humid),
                                          (u_grid,v_grid,vor_grid,div_grid,temp_grid,humid_grid))
                                           if length(M) > 0)
                                                    
            RingGrids.Matrix!(MGs...; quadrant_rotation, matrix_quadrant)

        else                    # or interpolate onto a full grid
            :u in output_vars       && RingGrids.interpolate!(output_Grid(u),    u_grid, interpolator)
            :v in output_vars       && RingGrids.interpolate!(output_Grid(v),    v_grid, interpolator)
            :vor in output_vars     && RingGrids.interpolate!(output_Grid(vor),  vor_grid, interpolator)
            :div in output_vars     && RingGrids.interpolate!(output_Grid(div),  div_grid, interpolator)
            :temp in output_vars    && RingGrids.interpolate!(output_Grid(temp), temp_grid, interpolator)
            :humid in output_vars   && RingGrids.interpolate!(output_Grid(humid),humid_grid, interpolator)
        end

        # UNSCALE THE SCALED VARIABLES
        unscale!(vor,diagn.scale[]) # was vor*radius, back to vor
        unscale!(div,diagn.scale[]) # same
        temp .-= 273.15             # convert to ˚C

        # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
        :u in output_vars     && round!(u,    keepbits.u)
        :v in output_vars     && round!(v,    keepbits.v)
        :vor in output_vars   && round!(vor,  keepbits.vor)
        :div in output_vars   && round!(div,  keepbits.div)
        :temp in output_vars  && round!(temp, keepbits.temp)
        :humid in output_vars && round!(humid,keepbits.humid)

        # WRITE VARIABLES TO FILE, APPEND IN TIME DIMENSION
        :u in output_vars     && (netcdf_file["u"][:,:,k,i] = u)
        :v in output_vars     && (netcdf_file["v"][:,:,k,i] = v)
        :vor in output_vars   && (netcdf_file["vor"][:,:,k,i] = vor)
        :div in output_vars   && (netcdf_file["div"][:,:,k,i] = div)
        :temp in output_vars  && (netcdf_file["temp"][:,:,k,i] = temp)
        :humid in output_vars && (netcdf_file["humid"][:,:,k,i] = humid)
    end

    # surface pressure, i.e. interface displacement η
    (; pres_grid, precip_large_scale, precip_convection, cloud_top ) = diagn.surface

    if output.as_matrix
        if :pres in output_vars || :precip_cond in output_vars || :precip_conv in output_vars || :cloud in output_vars
            MGs = ((M,G) for (M,G) in zip((pres,precip_cond,precip_conv,cloud),
                                            (pres_grid, precip_large_scale, precip_convection, cloud_top)))

            RingGrids.Matrix!(MGs...; quadrant_rotation, matrix_quadrant)
        end
    else
        :pres in output_vars && RingGrids.interpolate!(output_Grid(pres),pres_grid,interpolator)
        :precip in output_vars && RingGrids.interpolate!(output_Grid(precip_cond),precip_large_scale,interpolator)
        :precip in output_vars && RingGrids.interpolate!(output_Grid(precip_conv),precip_convection,interpolator)
        :precip in output_vars && RingGrids.interpolate!(output_Grid(precip_conv),precip_convection,interpolator)
        :cloud in output_vars && RingGrids.interpolate!(output_Grid(cloud),cloud_top,interpolator)
    end

    # after output set precip accumulators back to zero
    precip_large_scale .= 0
    precip_convection .= 0

    # convert from [m] to [mm] within output time step (e.g. 6hours)
    precip_cond *= 1000
    precip_conv *= 1000

    if Model <: PrimitiveEquation
        @. pres = exp(pres)/100     # convert from log(pₛ) to surface pressure pₛ [hPa]
    end

    :pres in output_vars && round!(pres,keepbits.pres)
    :precip in output_vars && round!(precip_cond,keepbits.precip_cond)
    :precip in output_vars && round!(precip_conv,keepbits.precip_conv)
    :cloud in output_vars && round!(cloud,keepbits.cloud)
    
    :pres in output_vars && (netcdf_file["pres"][:,:,i] = pres)
    :precip in output_vars && (netcdf_file["precip_cond"][:,:,i] = precip_cond)
    :precip in output_vars && (netcdf_file["precip_conv"][:,:,i] = precip_conv)
    :cloud in output_vars && (netcdf_file["cloud_top"][:,:,i] = cloud)

    return nothing
end

Base.close(output::OutputWriter) = NCDatasets.close(output.netcdf_file)
Base.close(::Nothing) = nothing

"""
$(TYPEDSIGNATURES)
A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions. The prognostic variables
are bitrounded for compression and the 2nd leapfrog time step is discarded.
Variables in restart file are unscaled."""
function write_restart_file(progn::PrognosticVariables{T},
                            output::OutputWriter) where T
    
    (; run_path, write_restart, keepbits ) = output
    output.output || return nothing         # exit immediately if no output and
    write_restart || return nothing         # exit immediately if no restart file desired
    
    # COMPRESSION OF RESTART FILE
    for layer in progn.layers

        # copy over leapfrog 2 to 1
        copyto!(layer.timesteps[1].vor,layer.timesteps[2].vor)
        copyto!(layer.timesteps[1].div,layer.timesteps[2].div)
        copyto!(layer.timesteps[1].temp,layer.timesteps[2].temp)
        copyto!(layer.timesteps[1].humid,layer.timesteps[2].humid)

        # bitround 1st leapfrog step to output precision
        if T <: Base.IEEEFloat  # currently not defined for other formats...
            round!(layer.timesteps[1].vor,keepbits.vor)
            round!(layer.timesteps[1].div,keepbits.div)
            round!(layer.timesteps[1].temp,keepbits.temp)
            round!(layer.timesteps[1].humid,keepbits.humid)
        end

        # remove 2nd leapfrog step by filling with zeros
        fill!(layer.timesteps[2].vor,0)
        fill!(layer.timesteps[2].div,0)
        fill!(layer.timesteps[2].temp,0)
        fill!(layer.timesteps[2].humid,0)
    end

    # same for surface pressure
    copyto!(progn.surface.timesteps[1].pres,progn.surface.timesteps[2].pres)
    T <: Base.IEEEFloat && round!(progn.surface.timesteps[1].pres,keepbits.pres)
    fill!(progn.surface.timesteps[2].pres,0)

    jldopen(joinpath(run_path,"restart.jld2"),"w"; compress=true) do f
        f["prognostic_variables"] = progn
        f["version"] = output.pkg_version
        f["description"] = "Restart file created for SpeedyWeather.jl"
    end
end

"""
$(TYPEDSIGNATURES)
Returns the full path of the output file after it was created.
"""
get_full_output_file_path(output::OutputWriter) = joinpath(output.run_path, output.filename)

"""
$(TYPEDSIGNATURES)
Loads a `var_name` trajectory of the model `M` that has been saved in a netCDF file during the time stepping.
"""
function load_trajectory(var_name::Union{Symbol, String}, model::ModelSetup) 
    @assert model.output.output "Output is turned off"
    return Array(NCDataset(get_full_output_file_path(model.output))[string(var_name)])
end
