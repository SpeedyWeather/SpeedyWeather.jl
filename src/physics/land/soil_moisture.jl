abstract type AbstractSoilMoisture <: AbstractParameterization end

export NoSoilMoisture
struct NoSoilMoisture <: AbstractSoilMoisture end
NoSoilMoisture(SG::SpectralGrid) = NoSoilMoisture()
initialize!(soil::NoSoilMoisture, model::PrimitiveEquation) = nothing

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    soil::NoSoilMoisture,
    model::PrimitiveEquation,
)
    return nothing
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    soil::NoSoilMoisture,
    model::PrimitiveEquation,
)
    return nothing
end

export SeasonalSoilMoisture
@kwdef struct SeasonalSoilMoisture{NF, Grid} <: AbstractSoilMoisture

    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    "number of soil layers"
    nlayers::Int

    # READ CLIMATOLOGY FROM FILE
    "[OPTION] path to the folder containing the soil moisture file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of soil moisture"
    file::String = "soil_moisture.nc"

    "[OPTION] variable name in netcdf file"
    varname_layer1::String = "swl1"
    varname_layer2::String = "swl2"

    "[OPTION] Grid the soil moisture file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "[OPTION] The missing value in the data respresenting ocean"
    missing_value::NF = NF(NaN)

    # to be filled from file
    "Monthly soil moisture volume fraction [1], interpolated onto Grid"
    monthly_soil_moisture::Grid = zeros(Grid, nlat_half, nlayers, 12)
end

# generator function
function SeasonalSoilMoisture(SG::SpectralGrid; kwargs...)
    (; NF, GridVariable4D, nlayers_soil, nlat_half) = SG
    return SeasonalSoilMoisture{NF, GridVariable4D}(; nlat_half, nlayers=nlayers_soil, kwargs...)
end

# don't both initializing for the dry model
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
    nlat_half = RingGrids.get_nlat_half(soil.file_Grid, nx*ny)

    # the soil moisture from file but wrapped into a grid
    NF = eltype(monthly_soil_moisture)
    soil_moisture_file = zeros(soil.file_Grid{NF}, nlat_half, soil.nlayers, nt)

    for l in 1:nt   # read out monthly to swap dimensions to horizontal - vertical - time
        soil_moisture_file[:, 1, l] .= vec(ncfile[soil.varname_layer1].var[:, :, l])
        soil_moisture_file[:, 2, l] .= vec(ncfile[soil.varname_layer2].var[:, :, l])
    end

    fill_value1 = NF(ncfile[soil.varname_layer1].attrib["_FillValue"])
    fill_value2 = NF(ncfile[soil.varname_layer2].attrib["_FillValue"])
    fill_value1 === fill_value2 || @warn "Fill values are different for the two soil layers, use only from layer 1"
    soil_moisture_file[soil_moisture_file .=== fill_value1] .= soil.missing_value      # === to include NaN
    
    @boundscheck grids_match(monthly_soil_moisture, soil_moisture_file, vertical_only=true) ||
        throw(DimensionMismatch(monthly_soil_moisture, soil_moisture_file))

    # create interpolator from grid in file to grid used in model
    interp = RingGrids.interpolator(Float32, monthly_soil_moisture, soil_moisture_file)
    interpolate!(monthly_soil_moisture, soil_moisture_file, interp)
    return nothing
end

function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    soil::SeasonalSoilMoisture,
    model::PrimitiveEquation,
)
    # initialize land temperature by "running" the step at the current time
    timestep!(progn, diagn, soil, model)
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
    weight = convert(NF, Dates.days(time-Dates.firstdayofmonth(time))/Dates.daysinmonth(time))

    (; monthly_soil_moisture) = soil
    (; soil_moisture) = progn.land

    for k in eachlayer(soil_moisture)
        for ij in eachgridpoint(soil_moisture)
            soil_moisture[ij, k] = (1-weight) * monthly_soil_moisture[ij, k, this_month] +
                                    weight  * monthly_soil_moisture[ij, k, next_month]
        end
    end

    return nothing
end

export LandBucketMoisture

@kwdef struct LandBucketMoisture{NF} <: AbstractSoilMoisture
    "[OPTION] Time scale of vertical diffusion [s]"
    time_scale::Second = Day(2)

    "[OPTION] Fraction of top layer runoff that is put into layer below [1]"
    runoff_fraction::NF = 0.5

    "[OPTION] Initial soil moisture, volume fraction [1]"
    initial_moisture::NF = 0
    
    "[OPTION] Apply land-sea mask to NaN ocean-only points?"
    mask::Bool = false

    "Field capacity per meter soil [m], top layer, f = γz, set by land.temperature"
    f₁::Base.RefValue{NF} = Ref(zero(NF))

    "Field capacity per meter soil [m], lower layer, f = γz, set by land.temperature"
    f₂::Base.RefValue{NF} = Ref(zero(NF))
end

LandBucketMoisture(SG::SpectralGrid; kwargs...) = LandBucketMoisture{SG.NF}(; kwargs...)
function initialize!(soil::LandBucketMoisture, model::PrimitiveEquation)
    (; nlayers_soil) = model.spectral_grid
    @assert nlayers_soil == 2 "LandBucketMoisture only works with 2 soil layers "*
        "but spectral_grid.nlayers_soil = $nlayers_soil given. Ignoring additional layers."

    # set the field capacity given layer thickness and 
    γ = model.land.thermodynamics.field_capacity
    z₁ = model.land.geometry.layer_thickness[1]
    z₂ = model.land.geometry.layer_thickness[2]
    soil.f₁[] = γ*z₁
    soil.f₂[] = γ*z₂

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
    set!(progn.land.soil_moisture, soil.initial_moisture)
    soil.mask && mask!(progn.land.soil_moisture, model.land_sea_mask, :ocean)
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

    Pconv = diagn.physics.precip_rate_convection    # precipitation in [m/s]
    Plsc = diagn.physics.precip_rate_large_scale
    E = diagn.physics.land.evaporative_flux         # [kg/s/m²], divide by density for [m/s]
    R = diagn.physics.land.river_runoff             # diagnosed [m/s]

    @boundscheck grids_match(soil_moisture, Pconv, Plsc, E, R, horizontal_only=true) || throw(DimensionMismatch(soil_moisture, Pconv))
    @boundscheck size(soil_moisture, 2) == 2 || throw(DimensionMismatch)
    f₁, f₂ = soil.f₁[], soil.f₂[]
    p = soil.runoff_fraction        # fraction of top layer runoff put into lower layer
    τ⁻¹ = inv(convert(eltype(soil_moisture), Second(soil.time_scale).value))
    f₁_f₂ = f₁/f₂
    Δt_f₁ = Δt/f₁

    @inbounds for ij in eachgridpoint(soil_moisture)
        if mask[ij] > 0                         # at least partially land
            # precipitation (convection + large-scale) minus evaporation
            # river runoff only diagnostic, i.e. R=0 here but drain excess water below
            # convert to [m/s] by dividing by density
            F = Pconv[ij] + Plsc[ij] - E[ij]/ρ    # - R[ij]

            # vertical diffusion term between layers
            D = τ⁻¹*(soil_moisture[ij, 1] - soil_moisture[ij, 2])

            # Equation in 8.5.2.2 of the MITgcm users guide (Land package)
            soil_moisture[ij, 1] += Δt_f₁*F - Δt*D
            soil_moisture[ij, 2] += Δt*f₁_f₂*D

            # river runoff
            W₁ = soil_moisture[ij, 1]
            δW₁ = W₁ - min(W₁, 1)               # excess moisture in top layer
            soil_moisture[ij, 1] -= δW₁         # remove excess from top layer
            soil_moisture[ij, 2] += p*δW₁*f₁_f₂ # add fraction to lower layer
            R[ij] += Δt*(1-p)*δW₁*f₁            # accumulate river runoff [m] of top layer

            # remove excess water from lower layer (this disappears)
            soil_moisture[ij, 2] = min(soil_moisture[ij, 2], 1)
        end
    end
end