"""
Abstract super type for ocean models, which control the sea surface temperature
and sea ice concentration as boundary conditions to a SpeedyWeather simulation.
A new ocean model has to be defined as

    CustomOceanModel <: AbstractOcean

and can have parameters like CustomOceanModel{T} and fields. They need to extend
the following functions

    function initialize!(ocean_model::CustomOceanModel, model::PrimitiveEquation)
        # your code here to initialize the ocean model itself
        # you can use other fields from model, e.g. model.geometry
    end

    function initialize!(   
        ocean::PrognosticVariablesOcean,
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        ocean_model::CustomOceanModel,
        model::PrimitiveEquation,
    )
        # your code here to initialize the prognostic variables for the ocean
        # namely, ocean.sea_surface_temperature, ocean.sea_ice_concentration, e.g.
        # ocean.sea_surface_temperature .= 300      # 300K everywhere
    end

    function ocean_timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        ocean_model::CustomOceanModel,
        model::PrimitiveEquation,
    )
        # your code here to change the progn.ocean.sea_surface_temperature and/or
        # progn.ocean.sea_ice_concentration on any timestep
    end

Temperatures in ocean.sea_surface_temperature have units of Kelvin,
or NaN for no ocean. Note that neither sea surface temperature, land-sea mask
or orography have to agree. It is possible to have an ocean on top of a mountain.
For an ocean grid-cell that is (partially) masked by the land-sea mask, its value will
be (fractionally) ignored in the calculation of surface fluxes (potentially leading
to a zero flux depending on land surface temperatures). For an ocean grid-cell that is NaN
but not masked by the land-sea mask, its value is always ignored.
"""
abstract type AbstractOcean end

function Base.show(io::IO, O::AbstractOcean)
    println(io, "$(typeof(O)) <: AbstractOcean")
    keys = propertynames(O)
    print_fields(io, O, keys)
end

# function barrier for all oceans
function initialize!(   ocean::PrognosticVariablesOcean,
                        progn::PrognosticVariables,
                        diagn::DiagnosticVariables,
                        model::PrimitiveEquation)
    initialize!(ocean, progn, diagn, model.ocean, model)
end

# function barrier for all oceans
function ocean_timestep!(   progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            model::PrimitiveEquation)
    ocean_timestep!(progn, diagn, model.ocean, model)
end


## SEASONAL OCEAN CLIMATOLOGY
export SeasonalOceanClimatology

"""
Seasonal ocean climatology that reads monthly sea surface temperature
fields from file, and interpolates them in time regularly
(default every 3 days) to be stored in the prognostic variables.
Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct SeasonalOceanClimatology{NF, Grid<:AbstractGrid{NF}} <: AbstractOcean

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    "[OPTION] Path to the folder containing the sea surface temperatures, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] Filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting land"
    missing_value::NF = NaN

    # to be filled from file
    "Monthly sea surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Vector{Grid} = [zeros(Grid, nlat_half) for _ in 1:12]
end

# generator function
function SeasonalOceanClimatology(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    return SeasonalOceanClimatology{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(ocean::SeasonalOceanClimatology, model::PrimitiveEquation)
    load_monthly_climatology!(ocean.monthly_temperature, ocean)
end

function load_monthly_climatology!( 
    monthly::Vector{Grid},
    scheme;
    varname::String = scheme.varname
) where {Grid<:AbstractGrid}

    # LOAD NETCDF FILE
    if scheme.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", scheme.file)
    else
        path = joinpath(scheme.path, scheme.file)
    end
    ncfile = NCDataset(path)

    # create interpolator from grid in file to grid used in model
    nx, ny = ncfile.dim["lon"], ncfile.dim["lat"]
    npoints = nx*ny
    NF_file = typeof(ncfile[varname].attrib["_FillValue"])
    grid = scheme.file_Grid(zeros(NF_file, npoints))
    interp = RingGrids.interpolator(Float32, monthly[1], grid)

    # interpolate and store in monthly
    for month in 1:12
        this_month = ncfile[varname][:, :, month]
        ij = 0
        for j in 1:ny
            for i in 1:nx
                ij += 1
                x = this_month[i, j]
                grid[ij] = ismissing(x) ? scheme.missing_value : x
            end
        end
        interpolate!(monthly[month], grid, interp)
    end
    return nothing
end

function initialize!(   
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::SeasonalOceanClimatology,
    model::PrimitiveEquation,
)
    interpolate_monthly!(   ocean.sea_surface_temperature,
                            ocean_model.monthly_temperature,
                            progn.clock.time)
end
    
function ocean_timestep!(   progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            ocean_model::SeasonalOceanClimatology,
                            model::PrimitiveEquation)

    interpolate_monthly!(   progn.ocean.sea_surface_temperature,
                            ocean_model.monthly_temperature,
                            progn.clock.time)
    return nothing
end

function interpolate_monthly!(
    grid::Grid,
    monthly::Vector{Grid},
    time::DateTime,
) where Grid

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    weight = convert(eltype(grid), Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    @. grid = (1-weight) * monthly[this_month] +
                 weight  * monthly[next_month]

    return nothing
end

## CONSTANT OCEAN CLIMATOLOGY
export ConstantOceanClimatology

"""
Constant ocean climatology that reads monthly sea surface temperature
fields from file, and interpolates them only for the initial conditions
in time to be stored in the prognostic variables. It is therefore an
ocean from climatology but without a seasonal cycle that is constant in time.
To be created like

    ocean = SeasonalOceanClimatology(spectral_grid)

and the ocean time is set with initialize!(model, time=time).
Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct ConstantOceanClimatology <: AbstractOcean
    "[OPTION] path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting land"
    missing_value::Float64 = NaN
end

# generator function, just pass on kwargs
function ConstantOceanClimatology(SG::SpectralGrid; kwargs...)
    ConstantOceanClimatology(; kwargs...)
end

# nothing to initialize for model.ocean
initialize!(::ConstantOceanClimatology, ::PrimitiveEquation) = nothing

# initialize 
function initialize!(   
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::ConstantOceanClimatology,
    model::PrimitiveEquation,
)
    Grid = typeof(ocean.sea_surface_temperature)
    (; nlat_half) = ocean.sea_surface_temperature
    monthly = [zeros(Grid, nlat_half) for _ in 1:12]
    load_monthly_climatology!(monthly, ocean_model)
    interpolate_monthly!(ocean.sea_surface_temperature, monthly, progn.clock.time)
end

function ocean_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::ConstantOceanClimatology,
    model::PrimitiveEquation,
)
    return nothing
end

## CONSTANT OCEAN CLIMATOLOGY
export AquaPlanet

"""
AquaPlanet sea surface temperatures that are constant in time and longitude,
but vary in latitude following a coslat². To be created like

    ocean = AquaPlanet(spectral_grid, temp_equator=302, temp_poles=273)

Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct AquaPlanet{NF} <: AbstractOcean
    "[OPTION] Temperature on the Equator [K]"
    temp_equator::NF = 302

    "[OPTION] Temperature at the poles [K]"
    temp_poles::NF = 273
end

# generator function
AquaPlanet(SG::SpectralGrid; kwargs...) = AquaPlanet{SG.NF}(; kwargs...)

# nothing to initialize for AquaPlanet
initialize!(::AquaPlanet, ::PrimitiveEquation) = nothing

# initialize 
function initialize!(   
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::AquaPlanet,
    model::PrimitiveEquation,
)
    (; sea_surface_temperature) = ocean
    Te, Tp = ocean_model.temp_equator, ocean_model.temp_poles
    sst(λ, φ) = (Te - Tp)*cosd(φ)^2 + Tp
    set!(sea_surface_temperature, sst, model.geometry)
end

function ocean_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::AquaPlanet,
    model::PrimitiveEquation,
)
    return nothing
end

export SlabOcean

@kwdef struct SlabOcean{NF} <: AbstractOcean
    "[OPTION] Initial temperature on the equator [K]"
    temp_equator::NF = 302

    "[OPTION] Initial temperature at the poles [K]"
    temp_poles::NF = 273

    "[OPTION] Specific heat capacity of water [J/kg/K]"
    specific_heat_capacity::NF = 4184

    "[OPTION] Average mixed-layer depth [m]"
    mixed_layer_depth::NF = 10

    "[OPTION] Density of water [kg/m³]"
    density::NF = 1000

    "[DERIVED] Effective mixed-layer heat capacity [J/K/m²]"
    heat_capacity_mixed_layer::NF = specific_heat_capacity*mixed_layer_depth*density
end

# generator function
SlabOcean(SG::SpectralGrid; kwargs...) = SlabOcean{SG.NF}(; kwargs...)

# nothing to initialize for SlabOcean
initialize!(ocean_model::SlabOcean, model::PrimitiveEquation) = nothing

# initialize 
function initialize!(
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::SlabOcean,
    model::PrimitiveEquation,
)
    # initialize like AquaPlanet
    (; sea_surface_temperature) = ocean
    Te, Tp = ocean_model.temp_equator, ocean_model.temp_poles
    sst(λ, φ) = (Te - Tp)*cosd(φ)^2 + Tp
    set!(sea_surface_temperature, sst, model.geometry)
end

function ocean_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::SlabOcean,
    model::PrimitiveEquation,
)
    sst = progn.ocean.sea_surface_temperature
    flux = diagn.physics.surface_flux_heat
    C₀ = ocean_model.heat_capacity_mixed_layer
    Δt = model.time_stepping.Δt_sec
    (; mask) = model.land_sea_mask

    # flux is defined as positive downward
    # Euler forward step, mask land fluxes
    @. sst += (Δt/C₀)*flux*(1 - mask)

    return nothing
end