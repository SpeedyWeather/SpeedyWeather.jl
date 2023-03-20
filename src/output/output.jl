Base.@kwdef mutable struct Output{NF<:Union{Float32,Float64}} <: AbstractOutput
    # NF: output only in Float32/64

    output::Bool = false                    # output to netCDF?
    output_vars::Vector{Symbol}=[:none]     # vector of output variables as Symbols
    write_restart::Bool = false             # also write restart file if output==true?
    
    startdate::DateTime = DateTime(2000,1,1)
    timestep_counter::Int = 0               # time step counter
    output_counter::Int = 0                 # output step counter
    n_timesteps::Int = 0                    # number of time steps
    n_outputsteps::Int = 0                  # number of time steps with output
    output_every_n_steps::Int = 0           # output every n time steps
    
    run_id::String = "-1"                   # run identification number
    run_path::String = ""                   # output path plus run????/

    # the netcdf file to be written into
    netcdf_file::Union{NcFile,Nothing} = nothing

    # grid specifications
    output_matrix::Bool = false             # if true sort grid points into a matrix (interpolation-free)
                                            # full grid for output if output_matrix == false
    output_Grid::Type{<:AbstractFullGrid} = FullGaussianGrid   
    nlat_half::Int = 0                      # size of the input/output grid
    nlon::Int = 0                           # number of longitude/latitude points
    nlat::Int = 0
    
    # input grid for interpolation
    input_Grid::Type{<:AbstractGrid} = FullGaussianGrid
    interpolator::AbstractInterpolator = DEFAULT_INTERPOLATOR(input_Grid,0,0)
    
    # fields to output (only one layer, reuse over layers)
    u::Matrix{NF} = zeros(0,0)              # zonal velocity
    v::Matrix{NF} = zeros(0,0)              # meridional velocity
    vor::Matrix{NF} = zeros(0,0)            # relative vorticity
    div::Matrix{NF} = zeros(0,0)            # divergence
    pres::Matrix{NF} = zeros(0,0)           # pressure
    temp::Matrix{NF} = zeros(0,0)           # temperature
    humid::Matrix{NF} = zeros(0,0)          # humidity
end

Output() = Output{Float32}()

"""
    run_id = get_run_id(output, output_path)

Checks existing `run-????` folders in output path to determine a 4-digit `run_id` number. 
"""
function get_run_id(output, output_path)

    if output 
        # pull list of existing run???? folders via readdir
        pattern = r"run-\d\d\d\d"                # run-???? in regex
        runlist = filter(x->startswith(x,pattern),readdir(output_path))
        runlist = filter(x->endswith(  x,pattern),runlist)
        existing_runs = [parse(Int,id[5:end]) for id in runlist]

        # get the run id from existing folders
        if length(existing_runs) == 0           # if no runfolder exists yet
            run_id = 1                          # start with run0001
        else
            run_id = maximum(existing_runs)+1   # next run gets id +1
        end
        
        return @sprintf("%04d",run_id)
    else
        return "-1"
    end
end

"""
    run_id, run_path = get_run_id_path(P::Parameters)

Creates a new folder `run-*` with the `run_id`. Also returns the full path
`run_path` of that folder. Returns `-1, "no runpath"` in the case of no output.
"""
function get_run_id_path(P::Parameters)

    @unpack output, output_path, run_id = P

    if output
        run_path = joinpath(output_path,string("run-",run_id_string(run_id)))
        @assert !(string("run-",run_id) in readdir(output_path)) "Run folder already exists, choose another run_id."

        mkdir(run_path)             # actually create the folder
        return run_id, run_path
    else
        return run_id, "no runpath"
    end
end

run_id_string(run_id::Integer) = @sprintf("%04d",run_id)
run_id_string(run_id::String) = run_id

"""
    outputter = initialize_netcdf_output(   progn::PrognosticVariables,
                                            diagn::DiagnosticVariables,
                                            M::ModelSetup)

Creates a netcdf file on disk and the corresponding `netcdf_file` object preallocated with output variables
and dimensions. `write_netcdf_output!` then writes consecuitive time steps into this file.
"""
function initialize_netcdf_output(  diagn::DiagnosticVariables,
                                    M::ModelSetup)

    M.parameters.output || return Output()          # escape directly when no netcdf output

    # DEFINE NETCDF DIMENSIONS TIME
    @unpack startdate = M.parameters
    time_string = "seconds since $(Dates.format(startdate, "yyyy-mm-dd HH:MM:0.0"))"
    dim_time = NcDim("time",0,unlimited=true)
    var_time = NcVar("time",dim_time,t=Int64,atts=Dict("units"=>time_string,"long_name"=>"time"))

    # DEFINE NETCDF DIMENSIONS SPACE
    @unpack output_matrix, output_Grid, output_nlat_half, nlev = M.parameters

    # if specified (>0) use output resolution via output_nlat_half, otherwise nlat_half from dynamical core
    nlat_half = output_nlat_half > 0 ? output_nlat_half : M.geometry.nlat_half

    if output_matrix == false   # interpolate onto (possibly different) output grid
        lond = get_lond(output_Grid,nlat_half)
        latd = get_latd(output_Grid,nlat_half)
        nlon = length(lond)
        nlat = length(latd)
        lon_name, lon_units, lon_longname = "lon","degrees_east","longitude"
        lat_name, lat_units, lat_longname = "lat","degrees_north","latitude"

    else                        # output grid directly into a matrix (resort grid points, no interpolation)
        @unpack nlat_half = M.geometry      # don't use output_nlat_half as not supported for output_matrix
        nlon,nlat = RingGrids.matrix_size(M.geometry.Grid,nlat_half)    # size of the matrix output
        lond = collect(1:nlon)                                          # just enumerate grid points for lond, latd
        latd = collect(1:nlat)
        lon_name, lon_units, lon_longname = "i","1","horizontal index i"
        lat_name, lat_units, lat_longname = "j","1","horizontal index j"
    end
    
    σ = M.geometry.σ_levels_full
    dim_lon = NcDim(lon_name,nlon,values=lond,atts=Dict("units"=>lon_units,"long_name"=>lon_longname))
    dim_lat = NcDim(lat_name,nlat,values=latd,atts=Dict("units"=>lat_units,"long_name"=>lat_longname))
    dim_lev = NcDim("lev",nlev,values=σ,atts=Dict("units"=>"1","long_name"=>"sigma levels"))

    # VARIABLES, define every variable here that could be output
    @unpack output_NF, compression_level = M.parameters
    missing_value = convert(output_NF,M.parameters.missing_value)

    # given pres the right name, depending on ShallowWaterModel or PrimitiveEquationModel
    pres_name = M isa ShallowWaterModel ? "interface displacement" : "surface pressure"
    pres_unit = M isa ShallowWaterModel ? "m" : "hPa"

    all_ncvars = (        # define NamedTuple to identify the NcVars by name
    time = var_time,
    u = NcVar("u",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"zonal wind","units"=>"m/s","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    v = NcVar("v",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"meridional wind","units"=>"m/s","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    vor = NcVar("vor",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"relative vorticity","units"=>"1/s","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    pres = NcVar("pres",[dim_lon,dim_lat,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>pres_name,"units"=>pres_unit,"missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    div = NcVar("div",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"divergence","units"=>"1/s","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    temp = NcVar("temp",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"temperature","units"=>"degC","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    humid = NcVar("humid",[dim_lon,dim_lat,dim_lev,dim_time],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"specific humidity","units"=>"1","missing_value"=>missing_value,
            "_FillValue"=>missing_value)),
    orog = NcVar("orog",[dim_lon,dim_lat],t=output_NF,compress=compression_level,
            atts=Dict("long_name"=>"orography","units"=>"m","missing_value"=>missing_value,
            "_FillValue"=>missing_value))
    )

    # CREATE NETCDF FILE
    @unpack output_filename, output_vars, output_NF = M.parameters
    run_id,run_path = get_run_id_path(M.parameters)     # create output folder and get its id and path

    # vector of NcVars for output
    vars_out = [all_ncvars[key] for key in keys(all_ncvars) if key in vcat(:time,output_vars)]
    netcdf_file = NetCDF.create(joinpath(run_path,output_filename),vars_out,mode=NetCDF.NC_NETCDF4)
    
    # CREATE OUTPUT STRUCT
    @unpack write_restart, startdate = M.parameters
    @unpack n_timesteps, n_outputsteps, output_every_n_steps = M.constants

    u = :u in output_vars ?         fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    v = :v in output_vars ?         fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    vor = :vor in output_vars ?     fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    div = :div in output_vars ?     fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    pres = :pres in output_vars ?   fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    temp = :temp in output_vars ?   fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)
    humid = :humid in output_vars ? fill(missing_value,nlon,nlat) : zeros(output_NF,0,0)

    # CREATE OUTPUT INTERPOLATOR
    input_Grid = M.parameters.Grid                      # grid to interpolate from (the on used in dyn core)
    Interpolator = M.parameters.output_Interpolator     # type of interpolator
    npoints = get_npoints(output_Grid,nlat_half)        # number of grid points to interpolate onto
    input_nlat_half = M.geometry.nlat_half
    interpolator = output_matrix ? Interpolator(input_Grid,0,0) :
                    Interpolator(output_NF,input_Grid,input_nlat_half,npoints)
    
    # PRECOMPUTE LOCATION INDICES
    latds, londs = get_latdlonds(output_Grid,nlat_half)
    output_matrix || update_locator!(interpolator,latds,londs)

    outputter = Output( output=true; output_vars, write_restart, 
                        startdate, n_timesteps, n_outputsteps,
                        output_every_n_steps, run_id, run_path,
                        netcdf_file, 
                        output_matrix, output_Grid, nlat_half, nlon, nlat,
                        input_Grid, interpolator,
                        u, v, vor, div, pres, temp, humid)

    # WRITE INITIAL CONDITIONS TO FILE
    write_netcdf_variables!(outputter,diagn,M)
    write_netcdf_time!(outputter,startdate)

    return outputter
end

"""write_netcdf_output!(outputter::Output,                      # contains everything for netcdf_file output
                        time_sec::Int,                          # model time [s] for output
                        progn::PrognosticVariables,             # all prognostic variables
                        diagn::DiagnosticVariables,             # all diagnostic variables
                        M::ModelSetup)                          # all parameters

Writes the variables from `diagn` of time step `i` at time `time_sec` into `netcdf_file`. Simply escapes for no
netcdf output of if output shouldn't be written on this time step. Converts variables from `diagn` to float32
for output, truncates the mantissa for higher compression and applies lossless compression."""
function write_netcdf_output!(  outputter::Output,              # everything for netcdf output
                                time::DateTime,                 # model time for output
                                diagn::DiagnosticVariables,     # all diagnostic variables
                                model::ModelSetup)              # all parameters

    outputter.timestep_counter += 1                                 # increase counter
    @unpack output, output_every_n_steps, timestep_counter = outputter
    output || return nothing                                        # escape immediately for no netcdf output
    timestep_counter % output_every_n_steps == 0 || return nothing  # escape if output not written on this step

    # WRITE VARIABLES
    write_netcdf_variables!(outputter,diagn,model)
    write_netcdf_time!(outputter,time)
end

function write_netcdf_time!(outputter::Output,
                            time::DateTime)
    
    @unpack netcdf_file, startdate = outputter
    i = outputter.output_counter

    time_sec = [round(Int64,Dates.value(Dates.Second(time-startdate)))]
    NetCDF.putvar(netcdf_file,"time",time_sec,start=[i])    # write time [sec] of next output step
    NetCDF.sync(netcdf_file)                                # sync to flush variables to disc

    return nothing
end

function write_netcdf_variables!(   outputter::Output,
                                    diagn::DiagnosticVariables,
                                    model::ModelSetup)

    outputter.output_counter += 1                   # increase output step counter
    @unpack output_vars = outputter                 # Vector{Symbol} of variables to output
    i = outputter.output_counter

    @unpack u, v, vor, div, pres, temp, humid = outputter
    @unpack output_matrix, output_Grid = outputter
    @unpack interpolator = outputter

    # output to matrix options
    quadrant_rotation = model.parameters.output_quadrant_rotation
    matrix_quadrant = model.parameters.output_matrix_quadrant

    for (k,diagn_layer) in enumerate(diagn.layers)
        
        @unpack u_grid, v_grid, vor_grid, div_grid, temp_grid, humid_grid = diagn_layer.grid_variables
        
        if output_matrix    # resort gridded variables interpolation-free into a matrix

            # create (matrix,grid) tuples for simultaneous grid -> matrix conversion  
            MGs = ((M,G) for (M,G) in zip((u,v,vor,div,temp,humid),
                                          (u_grid,v_grid,vor_grid,div_grid,temp_grid,humid_grid))
                                           if length(M) > 0)
                                                    
            RingGrids.Matrix!(MGs...; quadrant_rotation, matrix_quadrant)

        else                # or interpolate onto a full grid
            :u in output_vars       && RingGrids.interpolate!(output_Grid(u),    u_grid, interpolator)
            :v in output_vars       && RingGrids.interpolate!(output_Grid(v),    v_grid, interpolator)
            :vor in output_vars     && RingGrids.interpolate!(output_Grid(vor),  vor_grid, interpolator)
            :div in output_vars     && RingGrids.interpolate!(output_Grid(div),  div_grid, interpolator)
            :temp in output_vars    && RingGrids.interpolate!(output_Grid(temp), temp_grid, interpolator)
            :humid in output_vars   && RingGrids.interpolate!(output_Grid(humid),humid_grid, interpolator)
        end

        # UNSCALE THE SCALED VARIABLES
        unscale!(vor,model)     # was vor*radius, back to vor
        unscale!(div,model)     # same
        temp .-= 273.15         # convert to ˚C

        # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
        @unpack keepbits = model.parameters
        for var in (u,v,vor,div,temp,humid)
            round!(var,keepbits)
        end

        # WRITE VARIABLES TO FILE, APPEND IN TIME DIMENSION
        @unpack netcdf_file = outputter
        :u in output_vars     && NetCDF.putvar(netcdf_file,"u",    u,    start=[1,1,k,i],count=[-1,-1,1,1])
        :v in output_vars     && NetCDF.putvar(netcdf_file,"v",    v,    start=[1,1,k,i],count=[-1,-1,1,1])
        :vor in output_vars   && NetCDF.putvar(netcdf_file,"vor",  vor,  start=[1,1,k,i],count=[-1,-1,1,1])
        :div in output_vars   && NetCDF.putvar(netcdf_file,"div",  div,  start=[1,1,k,i],count=[-1,-1,1,1])
        :temp in output_vars  && NetCDF.putvar(netcdf_file,"temp", temp, start=[1,1,k,i],count=[-1,-1,1,1])
        :humid in output_vars && NetCDF.putvar(netcdf_file,"humid",humid,start=[1,1,k,i],count=[-1,-1,1,1])
    end

    # surface pressure, i.e. interface displacement η
    if :pres in output_vars
        @unpack pres_grid = diagn.surface

        if output_matrix
            Matrix!(pres,diagn.surface.pres_grid; quadrant_rotation, matrix_quadrant)
        else
            RingGrids.interpolate!(output_Grid(pres),pres_grid,interpolator)
        end

        if model isa PrimitiveEquation
            @. pres = exp(pres)/100     # convert from log(pₛ) to surface pressure pₛ [hPa]
        else
            round!(pres,model.parameters.keepbits)
        end
        
        NetCDF.putvar(outputter.netcdf_file,"pres",pres,start=[1,1,i],count=[-1,-1,1])
    end

    return nothing
end

"""
    write_restart_file( time::DateTime,
                        progn::PrognosticVariables,
                        outputter::Output,
                        M::ModelSetup)

A restart file `restart.jld2` with the prognostic variables is written
to the output folder (or current path) that can be used to restart the model.
`restart.jld2` will then be used as initial conditions. The prognostic variables
are bitround for compression and the 2nd leapfrog time step is discarded.
While the dynamical core may work with scaled variables, the restart file
contains these variables unscaled."""
function write_restart_file(time::DateTime,
                            progn::PrognosticVariables,
                            outputter::Output,
                            model::ModelSetup)
    
    @unpack run_path, write_restart = outputter
    write_restart || return nothing         # exit immediately if no restart file desired
    
    # COMPRESSION OF RESTART FILE
    @unpack keepbits = model.parameters
    for layer in progn.layers

        # copy over leapfrog 2 to 1
        copyto!(layer.leapfrog[1].vor,layer.leapfrog[2].vor)
        copyto!(layer.leapfrog[1].div,layer.leapfrog[2].div)
        copyto!(layer.leapfrog[1].temp,layer.leapfrog[2].temp)
        copyto!(layer.leapfrog[1].humid,layer.leapfrog[2].humid)

        # bitround 1st leapfrog step to output precision
        round!(layer.leapfrog[1].vor,keepbits)
        round!(layer.leapfrog[1].div,keepbits)
        round!(layer.leapfrog[1].temp,keepbits)
        round!(layer.leapfrog[1].humid,keepbits)

        # remove 2nd leapfrog step by filling with zeros
        fill!(layer.leapfrog[2].vor,0)
        fill!(layer.leapfrog[2].div,0)
        fill!(layer.leapfrog[2].temp,0)
        fill!(layer.leapfrog[2].humid,0)
    end

    # same for surface pressure
    copyto!(progn.pres.leapfrog[1],progn.pres.leapfrog[2])
    round!(progn.pres.leapfrog[1],keepbits)
    fill!(progn.pres.leapfrog[2],0)

    jldopen(joinpath(run_path,"restart.jld2"),"w"; compress=true) do f
        f["prognostic_variables"] = progn
        f["time"] = time
        f["version"] = model.parameters.version
        f["description"] = "Restart file created for SpeedyWeather.jl"
    end
end

"""
    get_full_output_file_path(p::Parameters)

Returns the full path of the output file after it was created.
"""
get_full_output_file_path(P::Parameters) = joinpath(P.output_path, string("run-",run_id_string(P.run_id),"/"), P.output_filename)

"""
    load_trajectory(var_name::Union{Symbol, String}, M::ModelSetup) 

Loads a `var_name` trajectory of the model `M` that has been saved in a netCDF file during the time stepping.
"""
load_trajectory(var_name::Union{Symbol, String}, M::ModelSetup) = NetCDF.ncread(get_full_output_file_path(M.parameters), string(var_name))
