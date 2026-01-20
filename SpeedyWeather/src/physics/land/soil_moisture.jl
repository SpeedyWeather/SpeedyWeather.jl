abstract type AbstractSoilMoisture <: AbstractLandComponent end

export SeasonalSoilMoisture

"""SeasonalSoilMoisture model that prescribes soil moisture from a monthly climatology file.
The soil moisture is linearly interpolated between months based on the model time.
$(TYPEDFIELDS)"""
@kwdef struct SeasonalSoilMoisture{NF, GridVariable4D} <: AbstractSoilMoisture
    # READ CLIMATOLOGY FROM FILE
    "[OPTION] path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of soil moisture"
    file::String = "soil_moisture.nc"

    "[OPTION] variable name in netcdf file for layer 1"
    varname_layer1::String = "swl1"

    "[OPTION] variable name in netcdf file for layer 2"
    varname_layer2::String = "swl2"

    "[OPTION] Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly soil moisture volume fraction [1], interpolated onto Grid"
    monthly_soil_moisture::GridVariable4D
end

# TODO to adapt create a ManualSeasonalSoilMoisture component like AlbedoClimatology is adapted to ManualAlbedo
# Adapt.adapt_structure(to, soil::SeasonalSoilMoisture) = adapt(to, ManualSeasonalSoilMoisture(soil.monthly_soil_moisture))

# generator function
function SeasonalSoilMoisture(SG::SpectralGrid, geometry::LandGeometry; kwargs...)
    (; NF, GridVariable4D, grid) = SG
    (; nlayers) = geometry
    monthly_soil_moisture = zeros(GridVariable4D, grid, nlayers, 12)
    return SeasonalSoilMoisture{NF, GridVariable4D}(; monthly_soil_moisture, kwargs...)
end

function variables(::SeasonalSoilMoisture)
    return (
        PrognosticVariable(name = :soil_moisture, dims = Grid3D(), desc = "Soil moisture content (fraction of capacity)", units = "1", namespace = :land),
    )
end

# don't bother initializing for the dry model
initialize!(soil::SeasonalSoilMoisture, model::PrimitiveDry) = nothing

function initialize!(soil::SeasonalSoilMoisture, model::PrimitiveEquation)
    (; monthly_soil_moisture) = soil

    # LOAD NETCDF FILE
    if soil.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../../input_data", soil.file)
    else
        path = joinpath(soil.path, soil.file)
    end
    ncfile = NCDataset(path)

    # read out netCDF data
    nx, ny, nt = ncfile.dim["lon"], ncfile.dim["lat"], ncfile.dim["time"]
    nlat_half = RingGrids.get_nlat_half(soil.file_Grid, nx * ny)
    grid = soil.file_Grid(nlat_half)

    # the soil moisture from file but wrapped into a grid
    NF = eltype(monthly_soil_moisture)
    nlayers = size(monthly_soil_moisture, 2)
    soil_moisture_file = zeros(NF, grid, nlayers, nt)

    for l in 1:nt   # read out monthly to swap dimensions to horizontal - vertical - time
        soil_moisture_file[:, 1, l] .= vec(ncfile[soil.varname_layer1].var[:, :, l])
        soil_moisture_file[:, 2, l] .= vec(ncfile[soil.varname_layer2].var[:, :, l])
    end

    fill_value1 = NF(ncfile[soil.varname_layer1].attrib["_FillValue"])
    fill_value2 = NF(ncfile[soil.varname_layer2].attrib["_FillValue"])
    fill_value1 === fill_value2 || @warn "Fill values are different for the two soil layers, use only from layer 1"
    soil_moisture_file[soil_moisture_file .=== fill_value1] .= soil.missing_value      # === to include NaN
    soil_moisture_file = on_architecture(model.architecture, soil_moisture_file)

    @boundscheck fields_match(monthly_soil_moisture, soil_moisture_file, vertical_only = true) ||
        throw(DimensionMismatch(monthly_soil_moisture, soil_moisture_file))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(monthly_soil_moisture, soil_moisture_file, NF = Float32)
    interpolate!(monthly_soil_moisture, soil_moisture_file, interp)
    return nothing
end

function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::SeasonalSoilMoisture,
        model::PrimitiveEquation,
    )
    # initialize soil_moisture by "running" the step at the current time
    return timestep!(progn, diagn, soil, model)
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::SeasonalSoilMoisture,
        model::PrimitiveDry,
    )
    return nothing
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::SeasonalSoilMoisture,
        model::PrimitiveEquation,
    )
    (; time) = progn.clock

    this_month = Dates.month(time)
    next_month = (this_month % 12) + 1      # mod for dec 12 -> jan 1

    # linear interpolation weight between the two months
    # TODO check whether this shifts the climatology by 1/2 a month
    NF = eltype(progn.land.soil_moisture)
    weight = convert(NF, Dates.days(time - Dates.firstdayofmonth(time)) / Dates.daysinmonth(time))

    (; monthly_soil_moisture) = soil
    (; soil_moisture) = progn.land

    launch!(
        architecture(soil_moisture), RingGridWorkOrder, size(soil_moisture),
        seasonal_soil_moisture_kernel!,
        soil_moisture, monthly_soil_moisture, weight, this_month, next_month
    )

    return nothing
end

@kernel inbounds = true function seasonal_soil_moisture_kernel!(
        soil_moisture, monthly_soil_moisture, weight, this_month, next_month
    )
    ij, k = @index(Global, NTuple)

    soil_moisture[ij, k] = (1 - weight) * monthly_soil_moisture[ij, k, this_month] +
        weight * monthly_soil_moisture[ij, k, next_month]
end

export LandBucketMoisture

"""LandBucketMoisture model with two soil layers exchanging moisture via vertical diffusion.
Forced by precipitation, evaporation, surface condensation, snow melt and river runoff drainage.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct LandBucketMoisture{NF} <: AbstractSoilMoisture
    "[OPTION] Time scale of vertical diffusion [s]"
    time_scale::Second = Day(2)

    "[OPTION] Infiltration fraction, that is, fraction of top layer runoff that is put into layer below [1]"
    @param infiltration_fraction::NF = 0.25 (bounds = 0 .. 1,)

    "[OPTION] Apply land-sea mask to set ocean-only points?"
    mask::Bool = true

    "[OPTION] Initial soil moisture over ocean, volume fraction [1]"
    @param ocean_moisture::NF = 0 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure LandBucketMoisture
LandBucketMoisture(SG::SpectralGrid, geometry::LandGeometry; kwargs...) = LandBucketMoisture{SG.NF}(; kwargs...)
function initialize!(soil::LandBucketMoisture, model::PrimitiveEquation)
    nlayers = get_soil_layers(model)
    @assert nlayers == 2 "LandBucketMoisture only works with 2 soil layers " *
        "but geometry.nlayers = $nlayers given. Ignoring additional layers."

    return nothing
end

function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::LandBucketMoisture,
        model::PrimitiveDry,
    )
    return nothing
end

function initialize!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::LandBucketMoisture,
        model::PrimitiveEquation,
    )
    # create a seasonal model, initialize it and the variables
    seasonal_model = SeasonalSoilMoisture(model.spectral_grid, model.land.geometry)
    initialize!(seasonal_model, model)
    initialize!(progn, diagn, seasonal_model, model)
    # (seasonal model will be garbage collected hereafter)

    # set ocean "soil" moisture points (100% ocean only)
    masked_value = soil.ocean_moisture
    return if soil.mask
        # TODO: broadcasting over views of Fields of GPUArrays doesn't work
        sm = progn.land.soil_moisture.data
        sm[isnan.(sm)] .= masked_value
        mask!(progn.land.soil_moisture, model.land_sea_mask, :ocean; masked_value)
    end
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::LandBucketMoisture,
        model::PrimitiveDry,
    )
    return nothing
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        soil::LandBucketMoisture,
        model::PrimitiveEquation,
    )
    (; soil_moisture) = progn.land
    Δt = model.time_stepping.Δt_sec
    ρ = model.atmosphere.water_density
    (; mask) = model.land_sea_mask

    P = diagn.physics.rain_rate                     # precipitation (rain only) in [m/s]
    S = diagn.physics.land.snow_melt_rate           # [kg/m²/s] includes snow runoff leakage water
    H = diagn.physics.land.surface_humidity_flux    # [kg/m²/s], divide by density for [m/s], positive up
    R = diagn.physics.land.river_runoff             # diagnosed here, accumulated [m]

    @boundscheck fields_match(soil_moisture, P, S, H, R, horizontal_only = true) ||
        throw(DimensionMismatch(soil_moisture, P))
    @boundscheck size(soil_moisture, 2) >= 2 || throw(DimensionMismatch)

    # Water at field capacity [m], top and lower layer γ*z₁ and γ*z₂
    γ = model.land.thermodynamics.field_capacity
    f₁ = γ * model.land.geometry.layer_thickness[1]
    f₂ = γ * model.land.geometry.layer_thickness[2]

    p = soil.infiltration_fraction      # Infiltration fraction: fraction of top layer runoff put into lower layer
    τ⁻¹ = inv(convert(eltype(soil_moisture), Second(soil.time_scale).value))

    params = (; ρ, Δt, f₁, f₂, p, τ⁻¹)  # pack into NamedTuple for kernel

    return launch!(
        architecture(soil_moisture), LinearWorkOrder, (size(soil_moisture, 1),),
        land_bucket_soil_moisture_kernel!, soil_moisture, mask, P, S, H, R, params
    )
end

@kernel inbounds = true function land_bucket_soil_moisture_kernel!(
        soil_moisture, mask, P, S, H, R, params
    )

    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land
        (; ρ, Δt, f₁, f₂, p, τ⁻¹) = params
        # precipitation (rain only, convection + large-scale) minus evaporation (or condensation, = humidity flux)
        # river runoff via drain excess water below (that's just gone)
        # convert to [m/s] by dividing by density

        # Soil top sources and sinks
        # note: rain water can increase soil moisture regardless of snow cover
        # [kg/m²/s] -> [m/s] for humidity flux H and snow melt rate S
        F = P[ij] + (S[ij] - H[ij]) / ρ

        # vertical diffusion term between layers
        D = τ⁻¹ * (soil_moisture[ij, 1] - soil_moisture[ij, 2])

        # Equation in 8.5.2.2 of the MITgcm users guide (Land package)
        soil_moisture[ij, 1] += Δt / f₁ * F - Δt * D
        soil_moisture[ij, 2] += Δt * f₁ / f₂ * D

        # river runoff
        W₁ = soil_moisture[ij, 1]           # wrt to field capacity so maximum is 1
        δW₁ = W₁ - min(W₁, 1)               # excess moisture in top layer, cap at field capacity
        soil_moisture[ij, 1] -= δW₁         # remove excess from top layer
        soil_moisture[ij, 2] += p * δW₁ * f₁ / f₂ # add fraction to lower layer
        R[ij] += Δt * (1 - p) * δW₁ * f₁            # accumulate river runoff [m] of top layer

        # remove excess water from lower layer (this disappears)
        soil_moisture[ij, 2] = min(soil_moisture[ij, 2], 1)
    end
end

function variables(::LandBucketMoisture)
    return (
        # Prognostic variables
        PrognosticVariable(name = :soil_moisture, dims = Grid3D(), desc = "Soil moisture content (fraction of capacity)", units = "1", namespace = :land),
        # Diagnostic variables written to diagn.physics
        DiagnosticVariable(name = :river_runoff, dims = Grid2D(), desc = "River runoff from soil moisture", units = "m/s", namespace = :land),
        # Diagnostic variables read from diagn.physics
        DiagnosticVariable(name = :rain_rate, dims = Grid2D(), desc = "Convective precipitation rate", units = "m/s"),
        DiagnosticVariable(name = :surface_humidity_flux, dims = Grid2D(), desc = "Surface humidity flux", units = "kg/s/m²", namespace = :land),
        DiagnosticVariable(name = :snow_melt_rate, dims = Grid2D(), desc = "Snow melt rate + snow runoff", units = "m/s", namespace = :land),
    )
end
