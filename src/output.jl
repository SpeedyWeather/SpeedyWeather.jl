"""Checks output folders to determine a 4-digit run id number."""
function get_run_id_path(P::Parameters)

    @unpack output,outpath = P

    if output
        runlist = filter(x->startswith(x,"run"),readdir(outpath))
        existing_runs = [parse(Int,id[4:end]) for id in runlist]
        run_id = maximum(existing_runs)+1
        runpath = joinpath(outpath,@sprintf("run%04d",run_id))
        mkdir(runpath)

        return run_id,runpath
    else
        return 0,"no runpath"
    end
end


function output_initialise( Diag::DiagnosticVariables,   # output grid-space variables only
                            G::GeoSpectral,
                            P::Parameters)

    @unpack nlon,nlat,nlev = G.geometry
    @unpack lon,lat = G.geometry
    @unpack output_startdate, compression_level = P

    # INITIALISE FOLDER
    run_id,runpath = get_run_id_path(P)
    file_name = @sprintf("run%04d",run_id)
    
    # TIME
    time_string = "hours since $(Dates.format(output_startdate, "yyyy-mm-dd HH:MM:0.0"))"
    timedim = NcDim("time",0,unlimited=true)
    timevar = NcVar("time",timedim,t=Int32,atts=Dict("units"=>time_string))

    # SPACE
    londim  = NcDim("lon",nlon,values=lon)
    latdim  = NcDim("lat",nlat,values=lat)
    levdim  = NcDim("lev",nlev,values=collect(1:nlev))
    lonvar  = NcVar("lon",londim,t=Float32,atts=Dict("long_name"=>"longitude"))
    latvar  = NcVar("lat",latdim,t=Float32,atts=Dict("long_name"=>"latitude"))
    levvar  = NcVar("lev",levdim,t=Int32,atts=Dict("long_name"=>"model level"))

    # VARIABLES
    uvar        = NcVar("u",[londim,latdim,levdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"zonal wind","units"=>"m/s","missing_value"=>-999999f0))
    vvar        = NcVar("v",[londim,latdim,levdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"meridional wind","units"=>"m/s","missing_value"=>-999999f0))
    Tvar        = NcVar("T",[londim,latdim,levdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"temperature","units"=>"K","missing_value"=>-999999f0))
    humidvar    = NcVar("humid",[londim,latdim,levdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"specific humidity","units"=>"1","missing_value"=>-999999f0))
    geopotvar   = NcVar("geopot",[londim,latdim,levdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"geopotential height","units"=>"m","missing_value"=>-999999f0))
    logp0var    = NcVar("logp0",[londim,latdim,timedim],t=Float32,compress=compression_level,
                    atts=Dict("long_name"=>"surface pressure","units"=>"Pa","missing_value"=>-999999f0))

    # CREATE NETCDF FILE
    netcdffile = NetCDF.create("$file_name",
                    [timevar,lonvar,latvar,levvar,uvar,vvar,Tvar,humidvar,geopotvar,logp0var],mode=NetCDF.NC_NETCDF4)

    # WRITE INITIAL CONDITIONS TO FILE
    write_output!(0,0,netcdffile,Diag,P)

    return netcdffile
end

function write_output!( iout::Int,
                        t::Real,
                        ncfile::NcFile,
                        Diag::DiagnosticVariables,
                        P::Parameters)

    @unpack u,v,Tabs,humid,logp0 = Diag.gridvars
    @unpack keepbits = P

    # always convert to float32 for output
    u_output = Float32.(u)
    v_output = Float32.(v)
    T_output = Float32.(Tabs)
    humid_output = Float32.(humid)
    logp0 = Float32.(logp0)

    # ROUNDING FOR ROUND+LOSSLESS COMPRESSION
    for var in (u_output,v_output,T_output,humid_output)
        round!(var,keepbits)
    end

    # WRITE TO FILE
    NetCDF.putvar(ncfile,"u",u_output,start=[1,1,1,iout],count=[-1,-1,-1,1])
    NetCDF.putvar(ncfile,"v",v_output,start=[1,1,1,iout],count=[-1,-1,-1,1])
    NetCDF.putvar(ncfile,"T",T_output,start=[1,1,1,iout],count=[-1,-1,-1,1])
    NetCDF.putvar(ncfile,"humid",humid_output,start=[1,1,1,iout],count=[-1,-1,-1,1])
    NetCDF.putvar(ncfile,"geopot",geopot_output,start=[1,1,1,iout],count=[-1,-1,-1,1])
    NetCDF.putvar(ncfile,"logp0",logp0_output,start=[1,1,iout],count=[-1,-1,1])

    # WRITE TIME
    NetCDF.putvar(ncfile,"time",Int32[t],start=[iout])
end