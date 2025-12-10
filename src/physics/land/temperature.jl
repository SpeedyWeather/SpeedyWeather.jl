abstract type AbstractLandTemperature <: AbstractParameterization end

export SeasonalLandTemperature
@kwdef mutable struct SeasonalLandTemperature{NF, GridVariable3D} <: AbstractLandTemperature
    "[OPTION] path to the folder containing the land temperature file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of land surface temperatures"
    file::String = "land_surface_temperature.nc"

    "[OPTION] variable name in netcdf file"
    varname::String = "lst"

    "[OPTION] Grid the land surface temperature file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NaN

    "[OPTION] Apply land-sea mask to use fallback ocean temperature for ocean-only points?"
    mask::Bool = true

    "[OPTION] Fallback ocean temperature when mask=true [K]"
    ocean_temperature::NF = 285

    # to be filled from file
    "Monthly land surface temperatures [K], interpolated onto Grid"
    monthly_temperature::GridVariable3D
end

# generator function
function SeasonalLandTemperature(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable3D, grid) = SG
    monthly_temperature = zeros(GridVariable3D, grid, 12)  # 12 months
    return SeasonalLandTemperature{NF, GridVariable3D}(; monthly_temperature, kwargs...)
end

function variables(::SeasonalLandTemperature)
    return (
        PrognosticVariable(name=:soil_temperature, dims=Grid3D(), namespace=:land),
    )
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
    fill_value = ncfile[land.varname].attrib["_FillValue"]
    lst = land.file_Grid(ncfile[land.varname].var[:, :, :], input_as=Matrix)
    lst[lst .=== fill_value] .= land.missing_value      # === to include NaN
    lst = on_architecture(model.architecture, lst)
    
    @boundscheck fields_match(monthly_temperature, lst, vertical_only=true) ||
        throw(DimensionMismatch(monthly_temperature, lst))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(monthly_temperature, lst, NF=Float32)
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

    launch!(architecture(soil_temperature), RingGridWorkOrder, size(soil_temperature),
            seasonal_land_temperature_kernel!,
            soil_temperature, monthly_temperature, weight, this_month, next_month)

    return nothing
end

@kernel inbounds=true function seasonal_land_temperature_kernel!(
    soil_temperature, monthly_temperature, weight, this_month, next_month
)
    ij, k  = @index(Global, NTuple)
    
    soil_temperature[ij, k] = (1 - weight) * monthly_temperature[ij, this_month] + 
                            weight * monthly_temperature[ij, next_month]
end

## CONSTANT LAND CLIMATOLOGY
export ConstantLandTemperature
@kwdef mutable struct ConstantLandTemperature{NF} <: AbstractLandTemperature
    "[OPTION] Globally constant temperature"
    temperature::NF = 285

    "[OPTION] Apply land-sea mask to NaN ocean-only points?"
    mask::Bool = true
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
    set!(progn.land.soil_temperature, land.temperature)
    land.mask && mask!(progn.land.soil_temperature, model.land_sea_mask, :ocean)
end

# temperature is constant so do nothing during land timestep
timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, land::ConstantLandTemperature, args...) = nothing

function variables(::ConstantLandTemperature)
    return (
        PrognosticVariable(name=:soil_temperature, dims=Grid3D(), namespace=:land),
    )
end

export LandBucketTemperature

"""MITgcm's two-layer soil model (https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html). Fields assert
$(TYPEDFIELDS)"""
@kwdef mutable struct LandBucketTemperature{NF} <: AbstractLandTemperature
    "[OPTION] Apply land-sea mask to set ocean-only points?"
    mask::Bool = true
    
    "[OPTION] Initial soil temperature over ocean [K]"
    ocean_temperature::NF = 285
end

# generator function
LandBucketTemperature(SG::SpectralGrid; kwargs...) = LandBucketTemperature{SG.NF}(; kwargs...)
function initialize!(land::LandBucketTemperature, model::PrimitiveEquation)
    (; nlayers_soil) = model.spectral_grid
    @assert nlayers_soil == 2 "LandBucketTemperature only works with 2 soil layers "*
        "but spectral_grid.nlayers_soil = $nlayers_soil given. Ignoring additional layers."
    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::LandBucketTemperature,
    model::PrimitiveEquation,
)
    # create a seasonal model, initialize it and the variables
    seasonal_model = SeasonalLandTemperature(model.spectral_grid)
    initialize!(seasonal_model, model)
    initialize!(progn, diagn, seasonal_model, model)
    # (seasonal model will be garbage collected hereafter)

    # set ocean "land" temperature points (100% ocean only)
    if land.mask
        masked_value = land.ocean_temperature
        # TODO currently requries .data because of broadcasting issues
        lst = progn.land.soil_temperature.data
        lst[isnan.(lst)] .= masked_value
        mask!(progn.land.soil_temperature, model.land_sea_mask, :ocean; masked_value)
    end
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    land::LandBucketTemperature,
    model::PrimitiveEquation,
)
    (; soil_temperature, soil_moisture) = progn.land
    Lᵥ = model.atmosphere.latent_heat_condensation
    Lᵢ = model.atmosphere.latent_heat_sublimation
    # ρ = model.atmosphere.water_density
    Δt = model.time_stepping.Δt_sec
    
    (; mask) = model.land_sea_mask
    (; thermodynamics, geometry) = model.land

    # Sum up flux F following Frierson et al. 2006, eq (1)
    # use separate land fluxes (not ocean)
    Rsd = diagn.physics.surface_shortwave_down          # before albedo reflection
    Rsu = diagn.physics.land.surface_shortwave_up       # only albedo reflection
    Rld = diagn.physics.surface_longwave_down           # all in [W/m²]
    Rlu = diagn.physics.land.surface_longwave_up
    S = diagn.physics.land.sensible_heat_flux

    # except these in [kg/s/m²]
    H = haskey(diagn.physics.land, :surface_humidity_flux) ? diagn.physics.land.surface_humidity_flux : nothing
    M = haskey(diagn.physics.land, :snow_melt_rate) ? diagn.physics.land.snow_melt_rate : nothing

    @boundscheck fields_match(soil_temperature, Rsd, Rsu, Rld, Rlu, S, horizontal_only=true) ||
        throw(DimensionMismatch(soil_temperature, Rs))
    @boundscheck size(soil_moisture, 2) == size(soil_temperature, 2) == 2 || throw(DimensionMismatch)
    
    λ = thermodynamics.heat_conductivity_dry_soil
    γ = thermodynamics.field_capacity
    Cw = thermodynamics.heat_capacity_water
    Cs = thermodynamics.heat_capacity_dry_soil
    z₁ = geometry.layer_thickness[1]
    z₂ = geometry.layer_thickness[2]
    Δ =  2λ/(z₁ + z₂)   # thermal diffusion operator [W/(m² K)]
    
    params = (; Lᵥ, Lᵢ, γ, Cw, Cs, z₁, z₂, Δ, Δt)

    launch!(architecture(soil_temperature), LinearWorkOrder, (size(soil_temperature, 1),),
        land_bucket_temperature_kernel!, soil_temperature, mask, soil_moisture, Rsd, Rsu, Rlu, Rld, S, H, M,
        params)

    return nothing
end

@kernel inbounds=true function land_bucket_temperature_kernel!(
    soil_temperature, mask, soil_moisture, Rsd, Rsu, Rlu, Rld, S, H, M, params)
    
    ij = @index(Global, Linear)

    if mask[ij] > 0                         # at least partially land
        (; Lᵥ, Lᵢ, γ, Cw, Cs, z₁, z₂, Δ, Δt) = params
        
        # Cooling from snow melt rate, in [W/m²] = [J/kg] * [kg/m²/s]
        Q_melt = isnothing(M) ? zero(Lᵢ) : Lᵢ * M[ij]

        # latent heat flux [W/m²], zero if H not available
        L = isnothing(H) ? zero(Lᵥ) : Lᵥ*H[ij]

        # total surface downward heat flux [W/m^2]
        F = Rsd[ij] - Rsu[ij] - Rlu[ij] + Rld[ij] - L - S[ij] - Q_melt

        # heat capacity of the (wet) soil layers 1 and 2 [J/(m³ K)]
        # ignore snow here
        C₁ = Cw * soil_moisture[ij, 1] * γ + Cs
        C₂ = Cw * soil_moisture[ij, 2] * γ + Cs

        # vertical diffusion term between layers
        D = Δ*(soil_temperature[ij, 1] - soil_temperature[ij, 2])

        # Equation in 8.5.2.2 of the MITgcm users guide (Land package)
        soil_temperature[ij, 1] += Δt/(z₁*C₁)*(F - D)
        soil_temperature[ij, 2] += Δt/(z₂*C₂)*D
    end
end

function variables(::LandBucketTemperature)
    return (
        # Prognostic variables
        PrognosticVariable(name=:soil_temperature, dims=Grid3D(), namespace=:land),
        PrognosticVariable(name=:soil_moisture, dims=Grid3D(), namespace=:land),
        # Diagnostic variables read from diagn.physics
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down", units="W/m²"),
        DiagnosticVariable(name=:surface_shortwave_up, dims=Grid2D(), desc="Surface shortwave radiation up", units="W/m²", namespace=:land),
        DiagnosticVariable(name=:surface_longwave_down, dims=Grid2D(), desc="Surface longwave radiation down", units="W/m²"),
        DiagnosticVariable(name=:surface_longwave_up, dims=Grid2D(), desc="Surface longwave radiation up", units="W/m²", namespace=:land),
        DiagnosticVariable(name=:sensible_heat_flux, dims=Grid2D(), desc="Sensible heat flux", units="W/m²", namespace=:land),
    )
end