@with_kw mutable struct Output{NF<:Union{Float32,Float64}}
    output::Bool = false                    # output to netCDF?
    output_vars::Vector{Symbol}=[:none]     # vector of
    write_restart::Bool = false             # also write restart file if output==true?
    i_out::Integer = 0                      # output step counter
    n_timesteps::Integer = 0                # number of time steps
    n_outputsteps::Integer = 0              # number of time steps with output
    output_every_n_steps::Integer = 0       # output every n time steps
    run_id::Integer = -1                    # run identification number
    run_path::String = ""                   # output path plus run????/

    # the netcdf file to be written into
    netcdf_file::Union{NcFile,Nothing} = nothing

    output_grid::Symbol = :full                     # or :matrix
    Grid::Type{<:AbstractGrid} = FullGaussianGrid   # full grid for output if output_grid == :full
    spectral_transform::SpectralTransform = SpectralTransform() # spectral transform for that grid

    # fields to output (only one layer, reuse over layers)
    u::Matrix{NF} = zeros(NF,0,0)           # zonal velocity
    v::Matrix{NF} = zeros(NF,0,0)           # meridional velocity
    vor::Matrix{NF} = zeros(NF,0,0)         # relative vorticity
    div::Matrix{NF} = zeros(NF,0,0)         # divergence
    pres::Matrix{NF} = zeros(NF,0,0)        # pressure
    temp::Matrix{NF} = zeros(NF,0,0)        # temperature
    humid::Matrix{NF} = zeros(NF,0,0)       # humidity
end

Output() = Output{Float32}()

mutable struct Feedback
    verbose::Bool                           # print feedback to REPL?
    output::Bool                            # store output (here: write to parameters and progress.txt?)

    # PROGRESS
    progress_meter::ProgressMeter.Progress  # struct containing everything progress related
    progress_txt::Union{IOStream,Nothing}   # txt is a Nothing in case of no output

    # NANS AND OTHER MODEL STATE FEEDBACK
    nans_detected::Bool                     # did NaNs occur in the simulation?
end

"""Initialises the progress txt file."""
function initialize_feedback(outputter::Output,M::ModelSetup)
    @unpack output, write_restart, n_timesteps = outputter
    @unpack run_id, run_path = outputter

    if output   # with netcdf output
        @unpack NF,n_days,trunc = M.parameters
        @unpack Grid,npoints = M.geometry

        # create progress.txt file in run????/
        progress_txt = open(joinpath(run_path,"progress.txt"),"w")
        s = "Starting SpeedyWeather.jl run $run_id on "*
                Dates.format(Dates.now(),Dates.RFC1123Format)
        write(progress_txt,s*"\n")      # and in file

        # add some information on resolution and number format
        write(progress_txt,"Integrating $(n_days) days at a spectral resolution of "*
                            "T$trunc on a $Grid with $npoints grid points.\n")
        write(progress_txt,"Number format is "*string(NF)*".\n")
        write(progress_txt,"All data will be stored in $run_path.\n")

        # also export parameters into run????/parameters.txt
        parameters_txt = open(joinpath(run_path,"parameters.txt"),"w")
        print(parameters_txt,M.parameters)
        close(parameters_txt)

    else        # no netcdf output
        progress_txt = nothing      # for no ouput, allocate dummies for Feedback struct
    end

    nans_detected = false           # currently not used

    # PROGRESSMETER
    @unpack verbose = M.parameters
    DT_IN_SEC[1] = M.constants.Δt_sec       # hack: redefine element in global constant dt_in_sec
                                            # used to pass on the time step to ProgressMeter.speedstring
    desc = "Weather is speedy$(output ? " run $run_id: " : ": ")"
                                            # show progress meter via `enabled` through verbose parameter
    progress_meter = ProgressMeter.Progress(n_timesteps,enabled=verbose,showspeed=true;desc)

    return Feedback(verbose,output,progress_meter,progress_txt,nans_detected)
end

"""Calls the progress meter and writes every 5% progress increase to txt."""
function progress!(F::Feedback)
    ProgressMeter.next!(F.progress_meter)           # update progress meter
    @unpack counter,n = F.progress_meter            # unpack counter after update

    # write progress to txt file too
    if (counter/n*100 % 1) > ((counter+1)/n*100 % 1)  
        percent = round(Int,(counter+1)/n*100)      # % of time steps completed
        if F.output && (percent % 5 == 0)           # write every 5% step in txt 
            write(F.progress_txt,"\n$percent%")
            flush(F.progress_txt)
        end
    end
end

"""Finalises the progress meter and the progress txt file."""
function progress_finish!(F::Feedback)
    ProgressMeter.finish!(F.progress_meter)
    
    if F.output     # write final progress to txt file
        time_elapsed = F.progress_meter.tlast - F.progress_meter.tinit
        s = "Integration done in $(readable_secs(time_elapsed))."
        write(F.progress_txt,"\n$s\n")
        flush(F.progress_txt)
    end
end

"""
    run_id, run_path = get_run_id_path(P::Parameters)

Checks existing `run????` folders in output path to determine a 4-digit `run_id` number
and creates a new folder `run????` with that `run_id`. Also returns the full path
`run_path` of that folder. Returns `0,"no runpath"` in the case of no output."""
function get_run_id_path(P::Parameters)

    @unpack output,output_path = P

    if output
        # pull list of existing run???? folders via readdir
        pattern = r"run\d\d\d\d"                # run???? in regex
        runlist = filter(x->startswith(x,pattern),readdir(output_path))
        runlist = filter(x->endswith(  x,pattern),runlist)
        existing_runs = [parse(Int,id[4:end]) for id in runlist]

        # get the run id from existing folders
        if length(existing_runs) == 0           # if no runfolder exists yet
            run_id = 1                          # start with run0001
        else
            run_id = maximum(existing_runs)+1   # next run gets id +1
        end

        run_path = joinpath(output_path,@sprintf("run%04d",run_id))
        mkdir(run_path)             # actually create the folder
        return run_id,run_path
    else
        return -1,"no runpath"
    end
end

"""
    outputter = initialize_netcdf_output(   progn::PrognosticVariables,
                                            diagn::DiagnosticVariables,
                                            M::ModelSetup)

Creates a netcdf file on disk and the corresponding `netcdf_file` object preallocated with output variables
and dimensions. `write_netcdf_output!` then writes consecuitive time steps into this file.
"""
function initialize_netcdf_output(  progn::PrognosticVariables,
                                    diagn::DiagnosticVariables,
                                    M::ModelSetup)

    M.parameters.output || return Output()          # escape directly when no netcdf output

    # DEFINE NETCDF DIMENSIONS TIME
    @unpack output_startdate = M.parameters
    time_string = "hours since $(Dates.format(output_startdate, "yyyy-mm-dd HH:MM:0.0"))"
    dim_time = NcDim("time",0,unlimited=true)
    var_time = NcVar("time",dim_time,t=Int32,atts=Dict("units"=>time_string,"long_name"=>"time"))

    # DEFINE NETCDF DIMENSIONS SPACE
    @unpack output_grid, nlev = M.parameters

    # Option :full, output via spectral transform on a full grid with lon,lat vectors
    # Option :matrix, output grid directly into a matrix (resort grid points, no interpolation)
    if output_grid == :full
        Grid = full_grid(M.geometry.Grid)           # use the full grid with same latitudes for output
        G = Geometry(M.parameters,Grid)         
        @unpack nlon,nlat = G                       # number of longitudes, latitudes
        @unpack lond,latd = G                       # lon, lat vectors in degree
        lon_name, lon_units, lon_longname = "lon","degrees_east","longitude"
        lat_name, lat_units, lat_longname = "lat","degrees_north","latitude"

    elseif output_grid == :matrix
        @unpack Grid, nresolution = M.geometry
        nlat,nlon = matrix_size(Grid,nresolution)   # size of the matrix output
        lond = collect(1:nlon)                      # use lond, latd, but just enumerate grid points
        latd = collect(1:nlat)
        lon_name, lon_units, lon_longname = "i","1","horizontal index i"
        lat_name, lat_units, lat_longname = "j","1","horizontal index j"
    end
    
    dim_lon = NcDim(lon_name,nlon,values=lond,atts=Dict("units"=>lon_units,"long_name"=>lon_longname))
    dim_lat = NcDim(lat_name,nlat,values=latd,atts=Dict("units"=>lat_units,"long_name"=>lat_longname))
    dim_lev = NcDim("lev",nlev,values=collect(Int32,1:nlev),atts=Dict("units"=>"1","long_name"=>"vertical model levels"))

    # VARIABLES, define every variable here that could be output
    @unpack output_NF, compression_level = M.parameters
    missing_value = convert(output_NF,M.parameters.missing_value)
    var_u   = NcVar("u",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"zonal wind","units"=>"m/s","missing_value"=>missing_value))
    var_v   = NcVar("v",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"meridional wind","units"=>"m/s","missing_value"=>missing_value))
    var_vor = NcVar("vor",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"relative vorticity","units"=>"1/s","missing_value"=>missing_value))
    var_pres    = NcVar("pres",[dim_lon,dim_lat,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"interface displacement","units"=>"m","missing_value"=>missing_value))
    var_div     = NcVar("div",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"divergence","units"=>"1/s","missing_value"=>missing_value))
    var_temp    = NcVar("temp",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"temperature","units"=>"K","missing_value"=>missing_value))
    var_humid   = NcVar("humid",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
        atts=Dict("long_name"=>"specific humidity","units"=>"1","missing_value"=>missing_value))

    # CREATE NETCDF FILE
    @unpack output_filename, output_vars, output_NF = M.parameters
    run_id,run_path = get_run_id_path(M.parameters)     # create output folder and get its id and path

    # vars_out = [var_time]               # vector of NcVars for output
    # for var in output_vars              # append var_* from above if specified in output_vars
    #     vari = Symbol(:var_,var)
    #     @eval push!($vars_out,$vari)
    # end
    vars_out = [var_time,var_vor]

    netcdf_file = NetCDF.create(joinpath(run_path,output_filename),vars_out,mode=NetCDF.NC_NETCDF4)
    
    # CREATE OUTPUT STRUCT
    @unpack write_restart = M.parameters
    @unpack n_timesteps, n_outputsteps, output_every_n_steps = M.constants
    spectral_transform = spectral_transform_for_full_grid(M.spectral_transform)

    u = :u in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    v = :v in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    vor = :vor in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    div = :div in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    pres = :pres in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    temp = :temp in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)
    humid = :humid in output_vars ? zeros(output_NF,nlon,nlat) : zeros(output_NF,0,0)

    outputter = Output( output=true; output_vars, write_restart, n_timesteps, n_outputsteps,
                        output_every_n_steps, run_id, run_path,
                        netcdf_file, output_grid, Grid, spectral_transform,
                        u, v, vor, div, pres, temp, humid)

    # WRITE INITIAL CONDITIONS TO FILE
    initial_time = 0                    # start at output_startdate
    lf = 1                              # output first leapfrog time step for initial conditions
    write_netcdf_variables!(outputter,progn,diagn,M,lf)
    write_netcdf_time!(outputter,initial_time)

    return outputter
end

"""write_netcdf_output!(netcdf_file::Union{NcFile,Nothing},     # netcdf file to output into
                        feedback::Feedback,                     # feedback struct to increment output counter
                        i::Int,                                 # time step index
                        time_sec::Int,                          # model time [s] for output
                        progn::PrognosticVariables,             # all prognostic variables
                        diagn::DiagnosticVariables,             # all diagnostic variables
                        M::ModelSetup)                          # all parameters

Writes the variables from `diagn` of time step `i` at time `time_sec` into `netcdf_file`. Simply escapes for no
netcdf output of if output shouldn't be written on this time step. Converts variables from `diagn` to float32
for output, truncates the mantissa for higher compression and applies lossless compression."""
function write_netcdf_output!(  outputter::Output,              # netcdf file to output into
                                feedback::Feedback,             # feedback struct to increment output counter
                                time_sec::Integer,              # model time [s] for output
                                progn::PrognosticVariables,     # all prognostic variables
                                diagn::DiagnosticVariables,     # all diagnostic variables
                                M::ModelSetup,                  # all parameters
                                lf::Integer=2)                  # leapfrog time step to output

    @unpack output, output_every_n_steps = outputter
    output || return nothing                                    # escape immediately for no netcdf output
    @unpack counter = feedback.progress_meter                   # feedback counts every time step
    counter % output_every_n_steps == 0 || return nothing       # escape if output not written on this step

    # WRITE VARIABLES
    write_netcdf_variables!(outputter,progn,diagn,M,lf)
    write_netcdf_time!(outputter,time_sec)
end

function write_netcdf_time!(outputter::Output,
                            time_sec::Integer)
    
    @unpack i_out, netcdf_file = outputter
    time_hrs = Int32[round(time_sec/3600)]                      # convert from seconds to hours
    NetCDF.putvar(netcdf_file,"time",time_hrs,start=[i_out])    # write time [hrs] of next output step
    NetCDF.sync(netcdf_file)                                    # sync to flush variables to disc

    return nothing
end

function write_netcdf_variables!(   outputter::Output,
                                    progn::PrognosticVariables,
                                    diagn::DiagnosticVariables,
                                    M::ModelSetup,
                                    lf::Integer=2)              # leapfrog time step to output

    outputter.i_out += 1                                        # increase output counter
    @unpack i_out, output_vars = outputter                      # Vector{Symbol} of variables to output

    @unpack u, v, vor, div, pres, temp, humid = outputter       # unpack all output variables
    @unpack output_grid, Grid = outputter                       # output_grid = :matrix or :full, Grid actual grid type
    S = outputter.spectral_transform                            # spectral transform for the (full) output grid

    for (k,(progn_layer,diagn_layer)) in enumerate(zip(progn.layers,diagn.layers))
        if output_grid == :matrix

            :u in output_vars       && Matrix!(u,     diagn_layer.grid_variables.U_grid)    # TODO somehow unscale coslat
            :v in output_vars       && Matrix!(v,     diagn_layer.grid_variables.V_grid)    # on the fly
            :vor in output_vars     && Matrix!(vor,   diagn_layer.grid_variables.vor_grid)
            :div in output_vars     && Matrix!(div,   diagn_layer.grid_variables.div_grid)
            :temp in output_vars    && Matrix!(temp,  diagn_layer.grid_variables.temp_grid)
            :humid in output_vars   && Matrix!(humid, diagn_layer.grid_variables.humid_grid)

        elseif output_grid == :full

            :u in output_vars       && gridded!(Grid(u),diagn_layer.dynamics_variables.u_coslat,S,unscale_coslat=true)
            :v in output_vars       && gridded!(Grid(v),diagn_layer.dynamics_variables.v_coslat,S,unscale_coslat=true)
            :vor in output_vars     && gridded!(Grid(vor),  progn_layer.leapfrog[lf].vor,S)
            :div in output_vars     && gridded!(Grid(div),  progn_layer.leapfrog[lf].div,S)
            :temp in output_vars    && gridded!(Grid(temp), progn_layer.leapfrog[lf].temp,S)
            :humid in output_vars   && gridded!(Grid(humid),progn_layer.leapfrog[lf].humid,S)
            
        end

        # UNSCALE SCALED VARIABLES
        @unpack radius_earth = M.geometry
        vor ./= M.geometry.radius_earth
        div ./= M.geometry.radius_earth

        # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
        @unpack keepbits = M.parameters
        for var in (u,v,vor,div,temp,humid)
            round!(var,keepbits)
        end

        # WRITE VARIABLES TO FILE, APPEND IN TIME DIMENSION
        @unpack netcdf_file = outputter
        :u in output_vars     && NetCDF.putvar(netcdf_file,"u",    u,    start=[1,1,k,i_out],count=[-1,-1,1,1])
        :v in output_vars     && NetCDF.putvar(netcdf_file,"v",    v,    start=[1,1,k,i_out],count=[-1,-1,1,1])
        :vor in output_vars   && NetCDF.putvar(netcdf_file,"vor",  vor,  start=[1,1,k,i_out],count=[-1,-1,1,1])
        :div in output_vars   && NetCDF.putvar(netcdf_file,"div",  div,  start=[1,1,k,i_out],count=[-1,-1,1,1])
        :temp in output_vars  && NetCDF.putvar(netcdf_file,"temp", temp, start=[1,1,k,i_out],count=[-1,-1,1,1])
        :humid in output_vars && NetCDF.putvar(netcdf_file,"humid",humid,start=[1,1,k,i_out],count=[-1,-1,1,1])
    end

    # surface pressure, i.e. interface displacement η
    if :pres in output_vars
        output_grid == :matrix && Matrix!(pres,diagn.pres_grid)
        output_grid == :full && gridded!(Grid(pres),progn.pres.leapfrog[lf],S)
        round!(pres,M.parameters.keepbits)
        NetCDF.putvar(netcdf_file,"pres",pres,start=[1,1,i_out],count=[-1,-1,1])
    end

    return nothing
end

"""
    write_restart_file( time_sec::Real,
                        progn::PrognosticVariables,
                        outputter::Output,
                        M::ModelSetup)

A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions. The prognostic variables
are bitround and the 2nd leapfrog time step is discarded."""
function write_restart_file(time_sec::Real,
                            progn::PrognosticVariables,
                            outputter::Output,
                            M::ModelSetup)
    
    @unpack run_path, write_restart = outputter
    write_restart || return nothing                 # exit immediately if no restart file desired

    # unscale variables
    @unpack radius_earth = M.geometry
    scale!(progn,:vor,1/radius_earth)
    scale!(progn,:div,1/radius_earth)

    # COMPRESSION OF RESTART FILE
    @unpack keepbits = M.parameters
    for layer in progn.layers
        # TODO copy from leapfrog 2 into 1 for storage

        # bitround 1st leapfrog step to output precision
        round!(layer.leapfrog[1].vor,keepbits)
        round!(layer.leapfrog[1].div,keepbits)
        round!(layer.leapfrog[1].temp,keepbits)
        round!(layer.leapfrog[1].humid,keepbits)

        # remove 2nd leapfrog step
        fill!(layer.leapfrog[2].vor,0)
        fill!(layer.leapfrog[2].div,0)
        fill!(layer.leapfrog[2].temp,0)
        fill!(layer.leapfrog[2].humid,0)
    end

    # same for surface pressure
    round!(progn.pres.leapfrog[1],keepbits)
    fill!(progn.pres.leapfrog[2],0)

    jldopen(joinpath(run_path,"restart.jld2"),"w"; compress=true) do f
        f["prognostic_variables"] = progn
        f["time"] = M.parameters.output_startdate + Dates.Second(time_sec)
        f["version"] = M.parameters.version
        f["description"] = "Restart file created for SpeedyWeather.jl"
    end
end