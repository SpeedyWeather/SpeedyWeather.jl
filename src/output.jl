"""
    run_id, run_path = get_run_id_path(P::Parameters)

Checks existing `run????` folders in output path to determine a 4-digit `run_id` number
and creates a new folder `run????` with that `run_id`. Also returns the full path
`run_path` of that folder. Returns `0,"no runpath"` in the case of no output."""
function get_run_id_path(P::Parameters)

    @unpack output,out_path = P

    if output
        # pull list of existing run???? folders via readdir
        pattern = r"run\d\d\d\d"                # run???? in regex
        runlist = filter(x->startswith(x,pattern),readdir(out_path))
        runlist = filter(x->endswith(  x,pattern),runlist)
        existing_runs = [parse(Int,id[4:end]) for id in runlist]

        # get the run id from existing folders
        if length(existing_runs) == 0           # if no runfolder exists yet
            run_id = 1                          # start with run0001
        else
            run_id = maximum(existing_runs)+1   # next run gets id +1
        end

        run_path = joinpath(out_path,@sprintf("run%04d",run_id))
        mkdir(run_path)             # actually create the folder
        return run_id,run_path
    else
        return 0,"no runpath"
    end
end

"""
    netcdf_file = initialize_netcdf_output( diagn::DiagnosticVariables, # output grid variables only
                                            feedback::Feedback,         # Feedback struct
                                            M::ModelSetup)              # ModelSetup struct

Creates a netcdf file on disk and the corresponding `netcdf_file` object preallocated with output variables
and dimensions. `write_netcdf_output!` then writes consecuitive time steps into this file.
"""
function initialize_netcdf_output(  diagn::DiagnosticVariables, # output grid variables only
                                    feedback::Feedback,         # Feedback struct
                                    M::ModelSetup)              # ModelSetup struct

    feedback.output || return nothing                   # escape directly when no netcdf output

    @unpack nlon,nlat,nlev = M.geospectral.geometry     # number of longitudes, latitudes, vertical levels
    @unpack lond,latd = M.geospectral.geometry          # lon, lat vectors in degree
    @unpack output_startdate, compression_level = M.parameters
    
    # DEFINE DIMENSIONS, TIME
    time_string = "hours since $(Dates.format(output_startdate, "yyyy-mm-dd HH:MM:0.0"))"
    dim_time = NcDim("time",0,unlimited=true)
    var_time = NcVar("time",dim_time,t=Int32,atts=Dict("units"=>time_string,"long_name"=>"time"))

    # AND SPACE
    dim_lon = NcDim("lon",nlon,values=lond,atts=Dict("units"=>"degrees_east","long_name"=>"longitude"))
    dim_lat = NcDim("lat",nlat,values=latd,atts=Dict("units"=>"degrees_north","long_name"=>"latitude"))
    dim_lev = NcDim("lev",nlev,values=collect(Int32,1:nlev),atts=Dict("units"=>"1","long_name"=>"vertical model levels"))

    # VARIABLES, u,v,vor are used for all ModelSetups
    var_u   = NcVar("u",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                atts=Dict("long_name"=>"zonal wind","units"=>"m/s","missing_value"=>-999999f0))
    var_v   = NcVar("v",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                atts=Dict("long_name"=>"meridional wind","units"=>"m/s","missing_value"=>-999999f0))
    var_vor = NcVar("vor",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                atts=Dict("long_name"=>"relative vorticity","units"=>"1/s","missing_value"=>-999999f0))

    # pressure is only used for ShallowWaterModel and PrimitiveEquationModel
    if typeof(M) <: ShallowWaterModel || typeof(M) <: PrimitiveEquationModel
        var_pres    = NcVar("pres",[dim_lon,dim_lat,dim_time],t=Float32,compress=compression_level,
                atts=Dict("long_name"=>"interface displacement","units"=>"m","missing_value"=>-999999f0))
        var_div     = NcVar("div",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                atts=Dict("long_name"=>"divergence","units"=>"1/s","missing_value"=>-999999f0))
    end

    # temperature and humidity only used for PrimitiveEquationModel
    if typeof(M) <: PrimitiveEquationModel
        var_temp    = NcVar("temp",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                        atts=Dict("long_name"=>"temperature","units"=>"K","missing_value"=>-999999f0))
        var_humid   = NcVar("humid",[dim_lon,dim_lat,dim_lev,dim_time],t=Float32,compress=compression_level,
                        atts=Dict("long_name"=>"specific humidity","units"=>"1","missing_value"=>-999999f0))
    end

    # CREATE NETCDF FILE
    @unpack run_id, run_path = feedback
    file_name = "output.nc"

    if typeof(M) <: BarotropicModel                 # output only u,v,vor
        netcdf_file = NetCDF.create(joinpath(run_path,file_name),[var_time,var_u,var_v,var_vor],
                        mode=NetCDF.NC_NETCDF4)
    elseif typeof(M) <: ShallowWaterModel           # output also divergence and pressure
        netcdf_file = NetCDF.create(joinpath(run_path,file_name),[var_time,var_u,var_v,var_vor,var_div,var_pres],
                        mode=NetCDF.NC_NETCDF4)
    elseif typeof(M) <: PrimitiveEquationModel      # output also temperature and humidity
        netcdf_file = NetCDF.create(joinpath(run_path,file_name),
                        [var_time,var_u,var_v,var_vor,var_div,var_pres,var_temp,var_humid],mode=NetCDF.NC_NETCDF4)
    end

    # WRITE INITIAL CONDITIONS TO FILE
    initial_time_sec = 0        # start at 0 hours after output_startdate
    write_netcdf_output!(netcdf_file,feedback,initial_time_sec,diagn,M)

    return netcdf_file
end

"""write_netcdf_output!(netcdf_file::Union{NcFile,Nothing},     # netcdf file to output into
                        feedback::Feedback,                     # feedback struct to increment output counter
                        i::Int,                                 # time step index
                        time_sec::Int,                          # model time [s] for output
                        diagn::DiagnosticVariables,             # all diagnostic variables
                        M::ModelSetup)                          # all parameters

Writes the variables from `diagn` of time step `i` at time `time_sec` into `netcdf_file`. Simply escapes for no
netcdf output of if output shouldn't be written on this time step. Converts variables from `diagn` to float32
for output, truncates the mantissa for higher compression and applies lossless compression."""
function write_netcdf_output!(  netcdf_file::Union{NcFile,Nothing},     # netcdf file to output into
                                feedback::Feedback,                     # feedback struct to increment output counter
                                time_sec::Int,                          # model time [s] for output
                                diagn::DiagnosticVariables,             # all diagnostic variables
                                M::ModelSetup)                          # all parameters

    @unpack counter = feedback.progress_meter
    isnothing(netcdf_file) && return nothing                            # escape immediately for no netcdf output
    counter % M.constants.output_every_n_steps == 0 || return nothing   # escape if output not written on this step

    feedback.i_out += 1                                         # increase counter
    @unpack i_out = feedback

    write_netcdf_variables!(i_out,netcdf_file,diagn,M)          # depending on ModelSetup M write variables to file

    # WRITE TIME
    time_hrs = Int32[round(time_sec/3600)]                      # convert from seconds to hours
    NetCDF.putvar(netcdf_file,"time",time_hrs,start=[i_out])    # write time [hrs] of next output step
    NetCDF.sync(netcdf_file)                                    # sync to flush variables to disc
end

function write_netcdf_variables!(   i_out::Integer,
                                    netcdf_file::NcFile,
                                    diagn::DiagnosticVariables,
                                    M::BarotropicModel)

    # CONVERT TO FLOAT32 FOR OUTPUT
    @unpack U_grid,V_grid,vor_grid = diagn.grid_variables
    u_output = convert.(Float32,U_grid)
    v_output = convert.(Float32,V_grid)
    vor_output = convert.(Float32,vor_grid)

    # UNSCALE SCALED VARIABLES
    unscale_coslat!(u_output,M.geospectral.geometry)
    unscale_coslat!(v_output,M.geospectral.geometry)
    vor_output ./= M.geospectral.geometry.radius_earth

    # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
    @unpack keepbits = M.parameters
    for var in (u_output,v_output,vor_output)
        round!(var,keepbits)
    end

    # WRITE VARIABLES TO FILE, APPEND IN TIME DIMENSION
    NetCDF.putvar(netcdf_file,"u",u_output,start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"v",v_output,start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"vor",vor_output,start=[1,1,1,i_out],count=[-1,-1,-1,1])
end

function write_netcdf_variables!(   i_out::Integer,
                                    netcdf_file::NcFile,
                                    diagn::DiagnosticVariables,
                                    M::ShallowWaterModel)

    # CONVERT TO FLOAT32 FOR OUTPUT
    @unpack U_grid,V_grid,vor_grid,div_grid,pres_grid = diagn.grid_variables
    u_output = convert.(Float32,U_grid)
    v_output = convert.(Float32,V_grid)
    vor_output = convert.(Float32,vor_grid)
    div_output = convert.(Float32,div_grid)
    pres_output = convert.(Float32,pres_grid)

    # UNSCALE SCALED VARIABLES
    unscale_coslat!(u_output,M.geospectral.geometry)
    unscale_coslat!(v_output,M.geospectral.geometry)
    vor_output ./= M.geospectral.geometry.radius_earth
    div_output ./= M.geospectral.geometry.radius_earth

    # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
    @unpack keepbits = M.parameters
    for var in (u_output,v_output,vor_output,div_output,pres_output)
        round!(var,keepbits)
    end

    # WRITE VARIABLES TO FILE, APPEND IN TIME DIMENSION
    NetCDF.putvar(netcdf_file,"u",  u_output,  start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"v",  v_output,  start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"vor",vor_output,start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"div",div_output,start=[1,1,1,i_out],count=[-1,-1,-1,1])
    NetCDF.putvar(netcdf_file,"pres",pres_output,start=[1,1,i_out],count=[-1,-1,1])
end