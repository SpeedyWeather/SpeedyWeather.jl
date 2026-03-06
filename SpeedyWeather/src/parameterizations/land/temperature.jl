abstract type AbstractLandTemperature <: AbstractLandComponent end

export SeasonalLandTemperature

"""SeasonalLandTemperature model that prescribes land surface temperature from a monthly climatology file.
The temperature is linearly interpolated between months based on the model time.
$(TYPEDFIELDS)"""
@kwdef struct SeasonalLandTemperature{NF, GridVariable3D} <: AbstractLandTemperature
    "[OPTION] filename of land surface temperatures"
    file::String = "land_surface_temperature.nc"

    "[OPTION] path to the folder containing the lst"
    path::String = joinpath("data", "boundary_conditions", file)

    "[OPTION] flag to check for lst in SpeedyWeatherAssets or locally"
    from_assets::Bool = true

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "[OPTION] variable name in netcdf file"
    varname::String = "lst"

    "[OPTION] Grid the land surface temperature file comes on"
    FieldType::Type{<:AbstractField} = FullGaussianField

    "[OPTION] Apply land-sea mask to use fallback ocean temperature for ocean-only points?"
    mask::Bool = true

    "[OPTION] Fallback ocean temperature when mask=true [K]"
    ocean_temperature::NF = 285

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::GridVariable3D
end

# TODO to adapt create a ManualSeasonalLandTemperature component like AlbedoClimatology is adapted to ManualAlbedo
# Adapt.adapt_structure(to, temp::SeasonalLandTemperature) = adapt(to, ManualSeasonalLandTemperature(temp.monthly_temperature))

# generator function
function SeasonalLandTemperature(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...)
    (; NF, GridVariable3D, grid) = SG
    monthly_temperature = zeros(GridVariable3D, grid, 12)  # 12 months
    return SeasonalLandTemperature{NF, GridVariable3D}(; monthly_temperature, kwargs...)
end

function variables(::SeasonalLandTemperature)
    return (
        PrognosticVariable(name = :soil_temperature, dims = Grid3D(), namespace = :land),
    )
end

function initialize!(land::SeasonalLandTemperature, model::PrimitiveEquation)
    (; monthly_temperature) = land

    # LOAD NETCDF FILE
    lst = get_asset(
        land.path;
        from_assets = land.from_assets,
        name = land.varname,
        ArrayType = land.FieldType,
        FileFormat = NCDataset,
        version = land.version
    )
    
    lst = on_architecture(model.architecture, lst)

    @boundscheck fields_match(monthly_temperature, lst, vertical_only = true) ||
        throw(DimensionMismatch(monthly_temperature, lst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(monthly_temperature, lst, NF = Float32)
    interpolate!(monthly_temperature, lst, interp)

    # mask ocean points to fallback ocean temperature
    # set ocean "land" temperature points (100% ocean only)
    masked_value = land.ocean_temperature
    if land.mask
        # Replace NaN values in soil temperature with a fallback ocean temperature
        # unpack via .data due to broadcasting issues
        mtd = monthly_temperature.data
        mtd[isnan.(mtd)] .= masked_value

        # but land-sea mask may not align so also set those 100% ocean points to
        # the same fallback ocean temperature
        mask!(monthly_temperature, model.land_sea_mask, :ocean; masked_value)
    end
    return nothing
end

function initialize!(
        vars::Variables,
        land::SeasonalLandTemperature,
        model::PrimitiveEquation,
    )
    # initialize land temperature by "running" the step at the current time
    return timestep!(vars, land, model)
end

function timestep!(
        vars::Variables,
        land::SeasonalLandTemperature,
        model::PrimitiveEquation,
    )
    (; time) = vars.prognostic.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    (; monthly_temperature) = land
    (; soil_temperature) = vars.prognostic.land
    NF = eltype(soil_temperature)
    weight = convert(NF, Dates.days(time - Dates.firstdayofmonth(time)) / Dates.daysinmonth(time))

    launch!(
        architecture(soil_temperature), RingGridWorkOrder, size(soil_temperature),
        seasonal_land_temperature_kernel!,
        soil_temperature, monthly_temperature, weight, this_month, next_month
    )

    return nothing
end

@kernel inbounds = true function seasonal_land_temperature_kernel!(
        soil_temperature, monthly_temperature, weight, this_month, next_month
    )
    ij, k = @index(Global, NTuple)

    soil_temperature[ij, k] = (1 - weight) * monthly_temperature[ij, this_month] +
        weight * monthly_temperature[ij, next_month]
end

## CONSTANT LAND CLIMATOLOGY
export ConstantLandTemperature
@parameterized @kwdef struct ConstantLandTemperature{NF} <: AbstractLandTemperature
    "[OPTION] Globally constant temperature"
    @param temperature::NF = 285 (bounds = Positive,)

    "[OPTION] Apply land-sea mask to NaN ocean-only points?"
    mask::Bool = true
end

# generator function
ConstantLandTemperature(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = ConstantLandTemperature{SG.NF}(; kwargs...)

initialize!(land::ConstantLandTemperature, model::PrimitiveEquation) = nothing
function initialize!(
        vars::Variables,
        land::ConstantLandTemperature,
        model::PrimitiveEquation,
    )
    set!(vars.prognostic.land.soil_temperature, land.temperature)
    land.mask && mask!(vars.prognostic.land.soil_temperature, model.land_sea_mask, :ocean)
    return nothing 
end

# temperature is constant so do nothing during land timestep
timestep!(vars::Variables, land::ConstantLandTemperature, args...) = nothing

function variables(::ConstantLandTemperature)
    return (
        PrognosticVariable(name = :soil_temperature, dims = Grid3D(), namespace = :land),
    )
end

export LandBucketTemperature

"""MITgcm's two-layer soil model (https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html).
Fields are $(TYPEDFIELDS)"""
@kwdef struct LandBucketTemperature{NF} <: AbstractLandTemperature
    "[OPTION] Apply land-sea mask to set ocean-only points?"
    mask::Bool = true

    "[OPTION] Initial soil temperature over ocean [K]"
    ocean_temperature::NF = 285
end

Adapt.@adapt_structure LandBucketTemperature

# generator function
LandBucketTemperature(SG::SpectralGrid, geometry::LandGeometryOrNothing = nothing; kwargs...) = LandBucketTemperature{SG.NF}(; kwargs...)
function initialize!(land::LandBucketTemperature, model::PrimitiveEquation)
    nlayers = get_nlayers(model.land)
    @assert nlayers == 2 "LandBucketTemperature only works with 2 soil layers " *
        "but geometry.nlayers = $nlayers given. Ignoring additional layers."
    return nothing
end

function initialize!(
        vars::Variables,
        land::LandBucketTemperature,
        model::PrimitiveEquation,
    )
    # create a seasonal model, initialize it and the variables
    seasonal_model = SeasonalLandTemperature(model.spectral_grid, model.land.geometry)
    initialize!(seasonal_model, model)
    initialize!(vars, seasonal_model, model)
    # (seasonal model will be garbage collected hereafter)

    # set ocean "land" temperature points (100% ocean only)
    if land.mask
        masked_value = land.ocean_temperature
        # TODO currently requries .data because of broadcasting issues
        lst = vars.prognostic.land.soil_temperature.data
        lst[isnan.(lst)] .= masked_value
        mask!(vars.prognostic.land.soil_temperature, model.land_sea_mask, :ocean; masked_value)
    end
    return nothing
end

function timestep!(
        vars::Variables,
        land::LandBucketTemperature,
        model::PrimitiveEquation,
    )
    (; soil_temperature, soil_moisture) = vars.prognostic.land
    Lᵥ = latent_heat_condensation(model.atmosphere)
    Lᵢ = latent_heat_sublimation(model.atmosphere)
    Δt = model.time_stepping.Δt_sec

    (; mask) = model.land_sea_mask
    (; thermodynamics, geometry) = model.land

    # Sum up flux F following Frierson et al. 2006, eq (1)
    # use separate land fluxes (not ocean)
    Rsd = vars.parameterizations.surface_shortwave_down         # before albedo reflection
    Rsu = vars.parameterizations.land.surface_shortwave_up      # only albedo reflection
    Rld = vars.parameterizations.surface_longwave_down          # all in [W/m²]
    Rlu = vars.parameterizations.land.surface_longwave_up
    S = vars.parameterizations.land.sensible_heat_flux

    # except these in [kg/s/m²]
    H = haskey(vars.parameterizations.land, :surface_humidity_flux) ? vars.parameterizations.land.surface_humidity_flux : nothing
    M = haskey(vars.parameterizations.land, :snow_melt_rate) ? vars.parameterizations.land.snow_melt_rate : nothing

    @boundscheck fields_match(soil_temperature, Rsd, Rsu, Rld, Rlu, S, horizontal_only = true) ||
        throw(DimensionMismatch(soil_temperature, Rs))
    @boundscheck size(soil_moisture, 2) == size(soil_temperature, 2) == 2 || throw(DimensionMismatch)

    λ = thermodynamics.heat_conductivity_dry_soil
    γ = thermodynamics.field_capacity
    Cw = thermodynamics.heat_capacity_water
    Cs = thermodynamics.heat_capacity_dry_soil
    z₁ = geometry.layer_thickness[1]
    z₂ = geometry.layer_thickness[2]
    Δ = 2λ / (z₁ + z₂)   # thermal diffusion operator [W/(m² K)]

    params = (; Lᵥ, Lᵢ, γ, Cw, Cs, z₁, z₂, Δ, Δt)

    launch!(
        architecture(soil_temperature), LinearWorkOrder, (size(soil_temperature, 1),),
        land_bucket_temperature_kernel!, soil_temperature, mask, soil_moisture, Rsd, Rsu, Rlu, Rld, S, H, M,
        params
    )

    return nothing
end

@kernel inbounds = true function land_bucket_temperature_kernel!(
        soil_temperature, mask, soil_moisture, Rsd, Rsu, Rlu, Rld, S, H, M, params
    )

    ij = @index(Global, Linear)

    if mask[ij] > 0                         # at least partially land
        (; Lᵥ, Lᵢ, γ, Cw, Cs, z₁, z₂, Δ, Δt) = params

        # Cooling from snow melt rate, in [W/m²] = [J/kg] * [kg/m²/s]
        Q_melt = isnothing(M) ? zero(Lᵢ) : Lᵢ * M[ij]

        # latent heat flux [W/m²], zero if H not available
        L = isnothing(H) ? zero(Lᵥ) : Lᵥ * H[ij]

        # total surface downward heat flux [W/m^2]
        F = Rsd[ij] - Rsu[ij] - Rlu[ij] + Rld[ij] - L - S[ij] - Q_melt

        # heat capacity of the (wet) soil layers 1 and 2 [J/(m³ K)]
        # ignore snow here
        C₁ = Cw * soil_moisture[ij, 1] * γ + Cs
        C₂ = Cw * soil_moisture[ij, 2] * γ + Cs

        # vertical diffusion term between layers
        D = Δ * (soil_temperature[ij, 1] - soil_temperature[ij, 2])

        # Equation in 8.5.2.2 of the MITgcm users guide (Land package)
        soil_temperature[ij, 1] += Δt / (z₁ * C₁) * (F - D)
        soil_temperature[ij, 2] += Δt / (z₂ * C₂) * D
    end
end

function variables(::LandBucketTemperature)
    return (
        PrognosticVariable(:soil_temperature, Grid3D(), desc = "Soil temperature", units = "K", namespace = :land),
        PrognosticVariable(:soil_moisture, Grid3D(), desc = "Soil moisture content (fraction of capacity)", units = "1", namespace = :land),
        ParameterizationVariable(:surface_shortwave_down, Grid2D(), desc = "Surface shortwave radiation down", units = "W/m²"),
        ParameterizationVariable(:surface_shortwave_up, Grid2D(), desc = "Surface shortwave radiation up", units = "W/m²", namespace = :land),
        ParameterizationVariable(:surface_longwave_down, Grid2D(), desc = "Surface longwave radiation down", units = "W/m²"),
        ParameterizationVariable(:surface_longwave_up, Grid2D(), desc = "Surface longwave radiation up", units = "W/m²", namespace = :land),
        ParameterizationVariable(:sensible_heat_flux, Grid2D(), desc = "Sensible heat flux", units = "W/m²", namespace = :land),
    )
end
