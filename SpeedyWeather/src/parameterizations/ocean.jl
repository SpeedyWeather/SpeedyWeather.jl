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

    function initialize!(vars::Variables, ocean_model::CustomOceanModel, model::PrimitiveEquation)
        # your code here to initialize the prognostic variables for the ocean
        # namely, vars.prognostic.ocean.sea_surface_temperature, e.g.
        # vars.prognostic.ocean.sea_surface_temperature .= 300      # 300K everywhere
    end

    function timestep!(vars::Variables, ocean_model::CustomOceanModel, model::PrimitiveEquation)
        # your code here to change the vars.prognostic.ocean.sea_surface_temperature
    end

Temperatures in ocean.sea_surface_temperature have units of Kelvin,
or NaN for no ocean. Note that neither sea surface temperature, land-sea mask
or orography have to agree. It is possible to have an ocean on top of a mountain.
For an ocean grid-cell that is (partially) masked by the land-sea mask, its value will
be (fractionally) ignored in the calculation of surface fluxes (potentially leading
to a zero flux depending on land surface temperatures). For an ocean grid-cell that is NaN
but not masked by the land-sea mask, its value is always ignored.
"""
abstract type AbstractOcean <: AbstractModelComponent end

# variable that AbstractOcean requires
function variables(::AbstractOcean)
    return (
        PrognosticVariable(:sea_surface_temperature, Grid2D(), namespace = :ocean, units = "K", desc = "Sea surface temperature"),
    )
end

# function barrier for all oceans
ocean_timestep!(vars::Variables, model::PrimitiveEquation) = timestep!(vars, model.ocean, model)
initialize!(vars::Variables, ocean::AbstractOcean, model) =
    @warn "`$(typeof(ocean))` does not have an `initialize!(::Variables, ::$(typeof(ocean)), model)` method " *
    "to set the initial conditions for the ocean prognostic variables."

export PrescribedOcean 

"""
Prescribed ocean that declares the necessary allocations for the sea surface temperature,
but all dynamics are expected to be externally set. Used for coupling to external ocean models.
"""
struct PrescribedOcean <: AbstractOcean end
PrescribedOcean(::SpectralGrid) = PrescribedOcean() # added constructor, just to be consistent with call signatures
initialize!(vars::Variables, ::PrescribedOcean, model) = nothing
timestep!(vars::Variables, ::PrescribedOcean, model) = nothing

export SeasonalOceanClimatology

"""
Seasonal ocean climatology that reads monthly sea surface temperature
fields from file, and interpolates them in time on every time step
and writes them to the prognostic variables.
Fields and options are
$(TYPEDFIELDS)"""
@kwdef struct SeasonalOceanClimatology{NF, Grid, GridVariable3D} <: AbstractOcean
    "Grid used for the model"
    grid::Grid

    "[OPTION] Filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] path to the folder containing the sst"
    path::String = joinpath("data", "boundary_conditions", file)

    "[OPTION] flag to check for sst in SpeedyWeatherAssets or locally"
    from_assets::Bool = true

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    FieldType::Type{<:AbstractField} = FullGaussianField

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

    sst = get_asset(
        ocean.path;
        from_assets = ocean.from_assets,
        name = ocean.varname,
        ArrayType = ocean.FieldType,
        FileFormat = NCDataset,
        version = ocean.version
    )

    # transfer to architecture of model if needed
    sst = on_architecture(model.architecture, sst)

    @boundscheck fields_match(monthly_temperature, sst, vertical_only = true) ||
        throw(DimensionMismatch(monthly_temperature, sst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(monthly_temperature, sst, NF = Float32)
    interpolate!(monthly_temperature, sst, interp)
    return nothing
end

# initial conditions for the ocean are just a "timestep" of the climatology
initialize!(vars::Variables, ocean::SeasonalOceanClimatology, model) = timestep!(vars, ocean, model)

function timestep!(
        vars::Variables,
        ocean::SeasonalOceanClimatology,
        model::PrimitiveEquation,
    )
    (; time) = vars.prognostic.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    (; monthly_temperature) = ocean
    (; sea_surface_temperature) = vars.prognostic.ocean
    NF = eltype(sea_surface_temperature)
    weight = convert(NF, Dates.days(time - Dates.firstdayofmonth(time)) / Dates.daysinmonth(year(time), Dates.month(time)))

    launch!(
        architecture(sea_surface_temperature), LinearWorkOrder, size(sea_surface_temperature),
        interpolate_monthly_climatology_kernel!,
        sea_surface_temperature, monthly_temperature, weight, this_month, next_month
    )
    return nothing
end

@kernel inbounds = true function interpolate_monthly_climatology_kernel!(
        var, monthly, weight, this_month, next_month
    )
    I = @index(Global, Cartesian)           # always launch over size(var)
    ijk = ndims(monthly) == 2 ? I[1] : I    # if monthly is 2D, ignore vertical index
    var[I] = (1 - weight) * monthly[ijk, this_month] + weight * monthly[ijk, next_month]
end

# CONSTANT OCEAN CLIMATOLOGY
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
    "[OPTION] filename of sea surface temperatures"
    file::String = "sea_surface_temperature.nc"

    "[OPTION] path to the folder containing the sst"
    path::String = joinpath("data", "boundary_conditions", file)

    "[OPTION] flag to check for sst in SpeedyWeatherAssets or locally"
    from_assets::Bool = true

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "[OPTION] Variable name in netcdf file"
    varname::String = "sst"

    "[OPTION] Grid the sea surface temperature file comes on"
    FieldType::Type{<:AbstractField} = FullGaussianField

    "[OPTION] The missing value in the data respresenting land"
    missing_value::Float64 = NaN
end

# generator
ConstantOceanClimatology(SG::SpectralGrid; kwargs...) = ConstantOceanClimatology(; kwargs...)

# nothing to initialize for model.ocean
initialize!(::ConstantOceanClimatology, ::PrimitiveEquation) = nothing

# initial conditions for the ocean are just a "timestep" of the climatology
function initialize!(vars::Variables, ocean_model::ConstantOceanClimatology, model)

    # create a seasonal model, initialize it and the variables
    (; path, file, varname, FieldType) = ocean_model
    (; NF, GridVariable3D, grid) = model.spectral_grid
    seasonal_model = SeasonalOceanClimatology{NF, typeof(grid), GridVariable3D}(;
        grid, path, file, varname, FieldType
    )
    initialize!(seasonal_model, model)

    # now set the initial conditions for the ocean prognostic variables with the seasonal model
    # (seasonal model will be garbage collected hereafter)
    initialize!(vars, seasonal_model, model)

    return nothing
end

# constant climatology does not change in time, so timestep! is a no-op
timestep!(vars::Variables, ocean_model::ConstantOceanClimatology, model) = nothing

# AQUAPLANET
export AquaPlanet

"""
AquaPlanet sea surface temperatures that are constant in time and longitude,
but vary in latitude following a coslat². To be created like

    ocean = AquaPlanet(spectral_grid, temp_equator=302, temp_poles=273)

Fields and options are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct AquaPlanet{NF} <: AbstractOcean
    "[OPTION] Temperature on the Equator [K]"
    @param temp_equator::NF = 302 (bounds = Positive,)

    "[OPTION] Temperature at the poles [K]"
    @param temp_poles::NF = 273 (bounds = Positive,)

    "[OPTION] Mask the sea surface temperature according to model.land_sea_mask?"
    mask::Bool = true
end

# generator function
AquaPlanet(SG::SpectralGrid; kwargs...) = AquaPlanet{SG.NF}(; kwargs...)

# nothing to initialize for AquaPlanet model itself
initialize!(::AquaPlanet, ::PrimitiveEquation) = nothing

# set initial conditions: cos²(lat) SST profile
function initialize!(vars::Variables, ocean_model::AquaPlanet, model::PrimitiveEquation)
    (; sea_surface_temperature) = vars.prognostic.ocean
    Te, Tp = ocean_model.temp_equator, ocean_model.temp_poles
    sst(λ, φ) = (Te - Tp) * cosd(φ)^2 + Tp
    set!(sea_surface_temperature, sst, model.geometry)
    ocean_model.mask && mask!(sea_surface_temperature, model.land_sea_mask, :land)
    return nothing
end

# SST is constant in time, so timestep! is a no-op
timestep!(vars::Variables, ocean_model::AquaPlanet, model::PrimitiveEquation) = nothing


export SlabOcean

@parameterized @kwdef mutable struct SlabOcean{NF} <: AbstractOcean
    "[OPTION] Specific heat capacity of water [J/kg/K]"
    specific_heat_capacity::NF = 4184

    "[OPTION] Average mixed-layer depth [m]"
    @param mixed_layer_depth::NF = 50 (bounds = Positive,)

    "[OPTION] Density of water [kg/m³]"
    density::NF = 1000

    "[OPTION] Mask initial sea surface temperature with land-sea mask?"
    mask::Bool = true

    "[OPTION] SST over land [K]"
    land_temperature::NF = 283

    "[DERIVED] Effective mixed-layer heat capacity [J/K/m²]"
    heat_capacity_mixed_layer::NF = specific_heat_capacity * mixed_layer_depth * density
end

# generator function
SlabOcean(SG::SpectralGrid; kwargs...) = SlabOcean{SG.NF}(; kwargs...)

function variables(::SlabOcean)
    return (
        PrognosticVariable(:sea_surface_temperature, Grid2D(), namespace = :ocean, desc = "Sea surface temperature", units = "K"),

        ParameterizationVariable(:surface_shortwave_down, Grid2D(), desc = "Surface shortwave radiation down", units = "W/m^2"),
        ParameterizationVariable(:surface_shortwave_up, Grid2D(), desc = "Surface shortwave radiation up over ocean", units = "W/m^2", namespace = :ocean),
        ParameterizationVariable(:surface_longwave_down, Grid2D(), desc = "Surface longwave radiation down", units = "W/m^2"),
        ParameterizationVariable(:surface_longwave_up, Grid2D(), desc = "Surface longwave radiation up over ocean", units = "W/m^2", namespace = :ocean),

        ParameterizationVariable(:surface_humidity_flux, Grid2D(), desc = "Surface humidity flux", units = "kg/s/m^2", namespace = :ocean),
        ParameterizationVariable(:surface_sensible_heat_flux, Grid2D(), desc = "Surface sensible heat flux", units = "kg/s/m^2", namespace = :ocean),
    )
end

# nothing to initialize for SlabOcean model itself
initialize!(ocean_model::SlabOcean, model::PrimitiveEquation) = nothing

# set initial conditions from seasonal climatology, then fill land points
function initialize!(vars::Variables, ocean_model::SlabOcean, model::PrimitiveEquation)
    # create a seasonal model, initialize it and the variables
    seasonal_model = SeasonalOceanClimatology(model.spectral_grid)
    initialize!(seasonal_model, model)
    initialize!(vars, seasonal_model, model)
    # (seasonal model will be garbage collected hereafter)

    # set land "sst" points (100% land only)
    if ocean_model.mask
        masked_value = ocean_model.land_temperature
        sst = vars.prognostic.ocean.sea_surface_temperature.data
        # TODO: broadcasting over views of Fields of GPUArrays doesn't work
        sst[isnan.(sst)] .= masked_value
        mask!(vars.prognostic.ocean.sea_surface_temperature, model.land_sea_mask, :land; masked_value)
    end
    return nothing
end

function timestep!(vars::Variables, ocean_model::SlabOcean, model::PrimitiveEquation)
    sst = vars.prognostic.ocean.sea_surface_temperature

    Lᵥ = latent_heat_condensation(model.atmosphere)
    C₀ = ocean_model.heat_capacity_mixed_layer
    Δt = model.time_stepping.Δt_sec
    Δt_C₀ = Δt / C₀

    (; mask) = model.land_sea_mask

    # Frierson et al. 2006, eq (1), all W/m² except humidity flux in kg/m²/s
    Rsd = vars.parameterizations.surface_shortwave_down         # before albedo
    Rsu = vars.parameterizations.ocean.surface_shortwave_up     # reflected from albedo
    Rld = vars.parameterizations.surface_longwave_down
    Rlu = vars.parameterizations.ocean.surface_longwave_up
    S = vars.parameterizations.ocean.sensible_heat_flux
    H = vars.parameterizations.ocean.surface_humidity_flux      # [kg/m²/s]

    params = (; Δt_C₀, Lᵥ)                              # pack into NamedTuple for kernel

    launch!(
        architecture(sst), LinearWorkOrder, size(sst), slab_ocean_kernel!,
        sst, mask, Rsd, Rsu, Rld, Rlu, H, S, params
    )
    return nothing
end

@kernel inbounds = true function slab_ocean_kernel!(sst, mask, Rsd, Rsu, Rld, Rlu, H, S, params)
    ij = @index(Global, Linear)         # every grid point ij
    if mask[ij] < 1                     # at least partially ocean
        (; Δt_C₀, Lᵥ) = params
        sst[ij] += Δt_C₀ * (Rsd[ij] - Rsu[ij] - Rlu[ij] + Rld[ij] - Lᵥ * H[ij] - S[ij])
    end
end
