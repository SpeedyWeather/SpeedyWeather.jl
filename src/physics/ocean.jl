"""
Abstract super type for ocean models, which control the sea surface temperature
and sea ice concentration as boundary conditions to a SpeedyWeather simulation.
A new ocean model has to be defined as

    CustomOceanModel <: AbstractOcean

and can have parameters like `CustomOceanModel{T}` and fields. They need to extend
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
    initialize!(ocean, progn, diagn, model.sea_ice, model)
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
@kwdef struct SeasonalOceanClimatology{NF, Grid, GridVariable3D} <: AbstractOcean

    "Grid used for the model"
    grid::Grid

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
    monthly_temperature::GridVariable3D = zeros(GridVariable3D, grid, 12)
end

# generator function
function SeasonalOceanClimatology(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable3D, grid) = SG
    return SeasonalOceanClimatology{NF, typeof(grid), GridVariable3D}(; grid, kwargs...)
end

function initialize!(ocean::SeasonalOceanClimatology, model::PrimitiveEquation)
    (; monthly_temperature) = ocean

    # LOAD NETCDF FILE
    if ocean.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", ocean.file)
    else
        path = joinpath(ocean.path, ocean.file)
    end
    ncfile = NCDataset(path)

    # create interpolator from grid in file to grid used in model
    fill_value = ncfile[ocean.varname].attrib["_FillValue"]
    sst = ocean.file_Grid(ncfile[ocean.varname].var[:, :, :], input_as=Matrix)
    sst[sst .=== fill_value] .= ocean.missing_value      # === to include NaN

    # transfer to architecture of model if needed 
    sst = on_architecture(model.architecture, sst)

    @boundscheck fields_match(monthly_temperature, sst, vertical_only=true) ||
        throw(DimensionMismatch(monthly_temperature, sst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(monthly_temperature, sst, NF=Float32)
    interpolate!(monthly_temperature, sst, interp)
    return nothing
end

function initialize!(
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::SeasonalOceanClimatology,
    model::PrimitiveEquation,
)
    ocean_timestep!(progn, diagn, ocean_model, model)
end

function ocean_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean::SeasonalOceanClimatology,
    model::PrimitiveEquation,
)
    (; time) = progn.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    (; monthly_temperature) = ocean
    (; sea_surface_temperature) = progn.ocean
    NF = eltype(sea_surface_temperature)
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    @inbounds for ij in eachgridpoint(sea_surface_temperature)
        sea_surface_temperature[ij] = (1-weight) * monthly_temperature[ij, this_month] +
                                          weight  * monthly_temperature[ij, next_month]
    end
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

and the ocean time is set with `initialize!(model, time=time)`.
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
    # create a seasonal model, initialize it and the variables
    (; path, file, varname, file_Grid, missing_value) = ocean_model
    (; NF, GridVariable3D, grid) = model.spectral_grid
    seasonal_model = SeasonalOceanClimatology{NF, typeof(grid), GridVariable3D}(;
                                grid, path, file, varname, file_Grid, missing_value)
    initialize!(seasonal_model, model)
    initialize!(ocean, progn, diagn, seasonal_model, model)
    # (seasonal model will be garbage collected hereafter)
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
    mask!(sea_surface_temperature, model.land_sea_mask, :land)
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

@kwdef mutable struct SlabOcean{NF, F} <: AbstractOcean
    "[OPTION] Initial temperature on the equator [K]"
    temp_equator::NF = 302

    "[OPTION] Initial temperature at the poles [K]"
    temp_poles::NF = 273

    "[OPTION] Specific heat capacity of water [J/kg/K]"
    specific_heat_capacity::NF = 4184

    "[OPTION] Average mixed-layer depth [m]"
    mixed_layer_depth::NF = 10

    "[OPTION] Sea ice insulation to reduce air-sea fluxes [0-1], 0->1 for no->full insulation"
    sea_ice_insulation::F

    "[OPTION] Density of water [kg/m³]"
    density::NF = 1000

    "[OPTION] Mask initial sea surface temperature with land-sea mask?"
    mask::Bool = false

    "[DERIVED] Effective mixed-layer heat capacity [J/K/m²]"
    heat_capacity_mixed_layer::NF = specific_heat_capacity*mixed_layer_depth*density
end

# generator function
function SlabOcean(
    SG::SpectralGrid;
    sea_ice_insulation = (x) -> x,  # default is linear reduction of air-sea fluxes with sea ice concentration
    kwargs...,
)
    return SlabOcean{SG.NF, typeof(sea_ice_insulation)}(; sea_ice_insulation, kwargs...)
end

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
    ocean_model.mask && mask!(sea_surface_temperature, model.land_sea_mask, :land)
end

function ocean_timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    ocean_model::SlabOcean,
    model::PrimitiveEquation,
)
    sst = progn.ocean.sea_surface_temperature
    ice = progn.ocean.sea_ice_concentration
    insulation = ocean_model.sea_ice_insulation

    Lᵥ = model.atmosphere.latent_heat_condensation
    C₀ = ocean_model.heat_capacity_mixed_layer
    Δt = model.time_stepping.Δt_sec
    Δt_C₀ = Δt / C₀

    (; mask) = model.land_sea_mask

    # Frierson et al. 2006, eq (1)
    Rsd = diagn.physics.surface_shortwave_down          # before albedo
    Rsu = diagn.physics.ocean.surface_shortwave_up      # reflected from albedo
    Rld = diagn.physics.surface_longwave_down
    Rlu = diagn.physics.ocean.surface_longwave_up
    Ev = diagn.physics.ocean.evaporative_flux
    S = diagn.physics.ocean.sensible_heat_flux

    launch!(architecture(sst), LinearWorkOrder, size(sst), slab_ocean_kernel!, sst, ice, mask, Rsd, Rsu, Rld, Rlu, Ev, S, insulation, Δt_C₀, Lᵥ)
    synchronize(architecture(sst))
end

@kernel inbounds=true function slab_ocean_kernel!(sst, ice, mask, Rsd, Rsu, Rld, Rlu, Ev, S, @Const(insulation), @Const(Δt_C₀), @Const(Lᵥ))
    ij = @index(Global, Linear)         # every grid point ij

    if mask[ij] < 1                     # at least partially ocean
        r = 1 - insulation(ice[ij])     # formulate as reduction of the net flux
        sst[ij] += Δt_C₀*r*(Rsd[ij] - Rsu[ij] - Rlu[ij] + Rld[ij] - Lᵥ*Ev[ij] - S[ij])
    end
end