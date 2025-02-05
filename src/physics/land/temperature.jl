abstract type AbstractLandTemperature <: AbstractParameterization end

export SeasonalLandTemperature
@kwdef struct SeasonalLandTemperature{NF, Grid} <: AbstractLand

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    # OPTIONS
    "[OPTION] path to the folder containing the land temperature file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of land surface temperatures"
    file::String = "land_surface_temperature.nc"

    "[OPTION] variable name in netcdf file"
    varname::String = "lst"

    "[OPTION] Grid the land surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::Grid = zeros(Grid, nlat_half, 12)
end

# generator function
function SeasonalLandTemperature(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable3D, nlat_half) = SG
    return SeasonalLandTemperature{NF, GridVariable3D}(; nlat_half, kwargs...)
end

function initialize!(land::SeasonalLandTemperature, model::PrimitiveEquation)
    (; monthly_temperature) = land

    # LOAD NETCDF FILE
    if land.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../../input_data", land.file)
    else
        path = joinpath(land.path, land.file)
    end
    ncfile = NCDataset(path)

    # read out netCDF data
    nx, ny, nt = ncfile.dim["lon"], ncfile.dim["lat"], ncfile.dim["time"]
    npoints = nx*ny
    fill_value = ncfile[land.varname].attrib["_FillValue"]
    lst = land.file_Grid(ncfile[land.varname].var[:, :, :], input_as=Matrix)
    lst[lst .=== fill_value] .= land.missing_value      # === to include NaN
    
    @boundscheck grids_match(monthly_temperature, lst, vertical_only=true) || throw(DimensionMismatch(monthly_temperature, lst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(Float32, monthly_temperature, lst)
    interpolate!(monthly_temperature, lst, interp)
    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::SeasonalLandTemperature,
    model::PrimitiveEquation,
)
    # initialize land temperature by "running" the step at the current time
    timestep!(progn, diagn, land, model)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::SeasonalLandTemperature,
    model::PrimitiveEquation,
)
    (; time) = progn.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    (; monthly_temperature) = land
    (; soil_temperature) = progn.land
    NF = eltype(soil_temperature)
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    k0 = 1   # first layer only
    for ij in eachgridpoint(soil_temperature)
        soil_temperature[ij, k0] = (1-weight) * monthly_temperature[ij, this_month] +
                                                weight  * monthly_temperature[ij, next_month]
    end

    # set other layers to the same temperature?
    # for k in eachgrid(soil_temperature)
    #     if k != k0
    #         soil_temperature[:, k] .= soil_temperature[:, k0]
    #     end
    # end

    return nothing
end

## CONSTANT LAND CLIMATOLOGY
export ConstantLandTemperature
@kwdef struct ConstantLandTemperature{NF} <: AbstractLand
    "[OPTION] Globally constant temperature"
    temperature::NF = 285
end

# generator function
ConstantLandTemperature(SG::SpectralGrid; kwargs...) = ConstantLandTemperature{SG.NF}(; kwargs...)

initialize!(land::ConstantLandTemperature, model::PrimitiveEquation) = nothing
function initialize!(   
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::ConstantLandTemperature,
    model::PrimitiveEquation,
)
    set!(progn.land.surface_temperature, land.temperature)
end

# temperature is constant so do nothing during land timestep
timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, land::ConstantLandTemperature, args...) = nothing

export LandBucketTemperature

"""MITgcm's two-layer soil model (https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html). Fields assert
$(TYPEDFIELDS)"""
@kwdef struct LandBucketTemperature{NF} <: AbstractLand
    "[OPTION] Top layer depth [m]"
    z₁::NF = 0.1

    "[OPTION] Second layer depth [m]"
    z₂::NF = 4.0

    "[OPTION] Thermal conductivity of the soil [W/(m K)]"
    λ::NF = 0.42

    "[OPTION] Field capacity per meter soil"
    γ::NF = 0.24

    "[OPTION] Heat capacity of water [J/(m³ K)]"
    Cw::NF = 4.2e6

    "[OPTION] Heat capacity of dry soil [J/(m³ K)]"
    Cs::NF = 1.13e6 

    "[OPTION] Initial soil temperature"
    initial_temperature::NF = 285

    "[OPTION] Apply land-sea mask to NaN ocean-only points?"
    mask::Bool = true
end

# generator function
LandBucketTemperature(SG::SpectralGrid; kwargs...) = LandBucketTemperature{SG.NF}(; kwargs...)
function initialize!(land::LandBucketTemperature, model::PrimitiveEquation)
    @assert model.spectral_grid.nlayers_soil < 2 "LandBucketTemperature only works with 2 soil layers"*
        "but spectral_grid.nlayers_soil = $(model.spectral_grid.nlayers_soil) given."
    @assert model.spectral_grid.nlayers_soil > 2 "LandBucketTemperature defined for 2 soil layers but"*
        "spectral_grid.nlayers_soil = $(model.spectral_grid.nlayers_soil) given. Ignoring additional layers."

    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::LandBucketTemperature,
    model::PrimitiveEquation,
)
    set!(progn.land.soil_temperature, land.initial_temperature)
    land.mask && mask!(progn.land.soil_temperature, model.land_sea_mask)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::LandBucketTemperature,
    model::PrimitiveEquation,
)
    (; soil_temperature, soil_moisture) = progn.land
    Lᵥ = model.atmosphere.latent_heat_condensation
    Δt = model.time_stepping.Δt_sec

    # Frierson et al. 2006, eq (1)
    Rs = diagn.physics.surface_shortwave_down
    Rld = diagn.physics.surface_longwave_down
    Rlu = diagn.physics.surface_longwave_up
    Ev = diagn.physics.evaporative_flux
    S = diagn.physics.sensible_heat_flux

    @boundscheck grids_match(soil_temperature, Rs, Rld, Rlu, Ev, S, horizontal_only=true) || throw(DimensionMismatch(soil_temperature, Rs))
    @boundscheck size(soil_moisture, 2) == size(soil_temperature, 2) == 2 || throw(DimensionMismatch)
    (; z₁, z₂, λ, γ, Cw, Cs) = land

    Δ =  2λ/(z₁ + z₂)   # thermal diffusion operator [W/(m² K)]

    for ij in eachgridpoint(soil_moisture, soil_temperature)

        # TODO if mask[ij] == at least partially land? only to skip ocean points?

        # total surface downward heat flux
        F = Rs[ij] - Rlu[ij] + Rld[ij] - Lᵥ*Ev[ij] - S[ij]

        # heat capacity of the soil layers 1 and 2 [J/(m³ K)]
        C₁ = Cw * soil_moisture[ij, 1] * γ + Cs
        C₂ = Cw * soil_moisture[ij, 2] * γ + Cs

        # vertical diffusion term between layers
        D = Δ*(soil_temperature[ij, 1] - soil_temperature[ij, 2])

        # Equation in 8.5.2.2 of the MITgcm users guide (Land package)
        soil_temperature[ij, 1] += Δt/(z₁*C₁)*(F - D)
        soil_temperature[ij, 2] += Δt/(z₂*C₂)*D
    end
end