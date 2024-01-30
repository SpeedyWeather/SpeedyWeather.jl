"""Coefficients to calculate the solar declination angle δ [radians] based on a simple
sine function, with Earth's axial tilt as amplitude, equinox as phase shift.
$(TYPEDFIELDS)"""
Base.@kwdef struct SinSolarDeclination{NF} <: AbstractSolarDeclination{NF}
    axial_tilt::NF = 23.44
    equinox::DateTime = DateTime(2000,3,20)
    length_of_year::Second = Day(365.25)
    length_of_day::Second = Hour(24)
end

"""Generator function pulling the number format NF from a SpectralGrid."""
SinSolarDeclination(SG::SpectralGrid;kwargs...) = SinSolarDeclination{SG.NF}(;kwargs...)

"""Generator function using the planet's orbital parameters to adapt the
solar declination calculation."""
function SinSolarDeclination(SG::SpectralGrid,P::AbstractPlanet)
    (;axial_tilt, equinox, length_of_year, length_of_day) = P
    SinSolarDeclination{SG.NF}(;axial_tilt, equinox, length_of_year, length_of_day)
end

"""
$(TYPEDSIGNATURES)
SinSolarDeclination functor, computing the solar declination angle of
angular fraction of year g [radians] using the coefficients of the
SinSolarDeclination struct."""
function (S::SinSolarDeclination)(g::NF) where NF
    axial_tilt = deg2rad(S.axial_tilt)
    equinox = S.length_of_day.value*Dates.dayofyear(S.equinox)/S.length_of_year.value
    return axial_tilt*sin(g-2*(π*convert(NF,equinox)))
end

"""Coefficients to calculate the solar declination angle δ from

    δ = 0.006918    - 0.399912*cos(g)  + 0.070257*sin(g)
                    - 0.006758*cos(2g) + 0.000907*sin(2g)
                    - 0.002697*cos(3g) + 0.001480*sin(3g)

with g the angular fraction of the year in radians. Following Spencer 1971,
Fourier series representation of the position of the sun. Search 2(5):172.
$(TYPEDFIELDS)"""
Base.@kwdef struct SolarDeclination{NF} <: AbstractSolarDeclination{NF}
    a::NF =   0.006918      # the offset +
    s1::NF =  0.070257      # s1*sin(g) +
    c1::NF = -0.399912      # c1*cos(g) +
    s2::NF =  0.000907      # s2*sin(2g) +
    c2::NF = -0.006758      # c2*cos(2g) +
    s3::NF =  0.001480      # s3*sin(3g) +
    c3::NF = -0.002697      # c3*cos(3g)
end

"""Generator function pulling the number format NF from a SpectralGrid."""
SolarDeclination(SG::SpectralGrid;kwargs...) = SolarDeclination{SG.NF}(;kwargs...)

"""
$(TYPEDSIGNATURES)
SolarDeclination functor, computing the solar declination angle of
angular fraction of year g [radians] using the coefficients of the
SolarDeclination struct."""
function (SD::SolarDeclination)(g)
    (;a,s1,s2,s3,c1,c2,c3) = SD
    sin1g,cos1g = sincos(g)
    sin2g,cos2g = sincos(2g)
    sin3g,cos3g = sincos(3g)
    return a + s1*sin1g + c1*cos1g + s2*sin2g + c2*cos2g + s3*sin3g + c3*cos3g
end

function Base.show(io::IO,L::AbstractSolarDeclination)
    println(io,"$(typeof(L)) <: AbstractSolarDeclination")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

"""Coefficients for the solar time correction (also called
Equation of time) which adjusts the solar hour to an oscillation
of sunrise/set by about +-16min throughout the year."""
Base.@kwdef struct SolarTimeCorrection{NF} <: AbstractSolarTimeCorrection{NF}
    a::NF =   0.004297      # the offset +
    s1::NF = -1.837877      # s1*sin(g) +
    c1::NF =  0.107029      # c1*cos(g) +
    s2::NF = -2.340475      # s2*sin(2g) +
    c2::NF = -0.837378      # c2*cos(2g)
end

"""
$(TYPEDSIGNATURES)
Functor that returns the time correction for a angular
fraction of the year g [radians], so that g=0 for Jan-01 and g=2π for Dec-31."""
function (STC::SolarTimeCorrection)(g)
    (;a,s1,s2,c1,c2) = STC
    sin1g,cos1g = sincos(g)
    sin2g,cos2g = sincos(2g)
    return deg2rad(a + s1*sin1g + c1*cos1g + s2*sin2g + c2*cos2g)
end

function Base.show(io::IO,L::AbstractSolarTimeCorrection)
    println(io,"$(typeof(L)) <: AbstractSolarTimeCorrection")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

"""
$(TYPEDSIGNATURES)
Chooses from SolarZenith (daily and seasonal cycle), SolarZenithDay (daily cycle only),
SolarZenithSeason (seasonal cycle only) and SolarZenithConstant (no daily or seasonal cycle)
given the parameters in model.planet."""
function WhichZenith(SG::SpectralGrid,P::AbstractPlanet;kwargs...)
    (;NF, Grid, nlat_half) = SG
    (;daily_cycle, seasonal_cycle, length_of_day, length_of_year) = P
    solar_declination = SinSolarDeclination(SG,P)

    if daily_cycle
        return SolarZenith{NF,Grid{NF}}(;
            nlat_half, length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...)

    else
        return SolarZenithSeason{NF,Grid{NF}}(;
            nlat_half, length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...)
    end
end

function cos_zenith!(time::DateTime,model::PrimitiveEquation)
    (;solar_zenith,geometry) = model
    cos_zenith!(solar_zenith,time,geometry)
end

function Base.show(io::IO,L::AbstractZenith)
    println(io,"$(typeof(L)) <: AbstractZenith")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

"""Solar zenith angle varying with daily and seasonal cycle.
$(TYPEDFIELDS)"""
Base.@kwdef struct SolarZenith{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractZenith{NF,Grid}
    # DIMENSIONS
    nlat_half::Int

    # OPTIONS
    length_of_day::Second = Hour(24)
    length_of_year::Second = Day(365.25)
    equation_of_time::Bool = true
    seasonal_cycle::Bool = true

    # COEFFICIENTS
    solar_declination::SinSolarDeclination{NF} = SinSolarDeclination{NF}()
    time_correction::SolarTimeCorrection{NF} = SolarTimeCorrection{NF}()

    # WORK ARRAYS
    initial_time::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
    cos_zenith::Grid = zeros(Grid,nlat_half)
end

function initialize!(
    S::AbstractZenith,
    initial_time::DateTime,
    model::ModelSetup
)
    S.initial_time[] = initial_time     # to fix the season if no seasonal cycle
    cos_zenith!(S,initial_time,model.geometry)
end

"""
$(TYPEDSIGNATURES)
Fraction of year as angle in radians [0...2π].
TODO: Takes length of day/year as argument, but calls to Dates.Time(), Dates.dayofyear()
currently have these hardcoded."""
function year_angle(::Type{T}, time::DateTime, length_of_day::Second, length_of_year::Second) where T
    year2rad = convert(T,2π/length_of_year.value)
    sec_of_day = Second(Dates.Time(time).instant).value
    return year2rad*(Dates.dayofyear(time)*length_of_day.value + sec_of_day)
end

"""
$(TYPEDSIGNATURES)
Fraction of day as angle in radians [0...2π].
TODO: Takes length of day as argument, but a call to Dates.Time()
currently have this hardcoded anyway."""
function solar_hour_angle(
    ::Type{T},
    time::DateTime,
    λ,                      # longitude in radians
    length_of_day::Second
) where T
    day2rad = convert(T,2π/length_of_day.value)
    noon_in_sec = length_of_day.value ÷ 2
    sec_of_day = Second(Dates.Time(time).instant).value
    return (sec_of_day - noon_in_sec)*day2rad + convert(T,λ)
end

function cos_zenith!(
    S::SolarZenith{NF},
    time::DateTime,
    geometry::Geometry
) where NF

    (;sinlat, coslat, lons) = geometry
    (;cos_zenith, length_of_day, length_of_year) = S

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? time : S.initial_time[]
    g = year_angle(NF, time_of_year, length_of_day, length_of_year)

    # time correction [radians] due to the equation of time (sunrise/set oscillation)
    tc = S.equation_of_time ? S.time_correction(g) : 0
    
    # solar hour angle at 0˚E (longtiude offset added later)
    λ = 0
    solar_hour_angle_0E = solar_hour_angle(NF, time, λ, length_of_day) + tc

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ,cosδ = sincos(δ)

    rings = eachring(cos_zenith)

    @inbounds for (j,ring) in enumerate(rings)                         
        sinδsinϕ = sinδ*sinlat[j]
        cosδcosϕ = cosδ*coslat[j]
        for ij in ring
            h = solar_hour_angle_0E + lons[ij]      # solar hour angle at longitude λ in radians
            cos_zenith[ij] = max(0,sinδsinϕ + cosδcosϕ*cos(h))
        end
    end
end

"""Solar zenith angle varying with seasonal cycle only.
$(TYPEDFIELDS)"""
Base.@kwdef struct SolarZenithSeason{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractZenith{NF,Grid}
    # DIMENSIONS
    nlat_half::Int

    # OPTIONS
    length_of_day::Second = Hour(24)
    length_of_year::Second = Day(365.25)
    seasonal_cycle::Bool = true

    # COEFFICIENTS
    solar_declination::SinSolarDeclination{NF} = SinSolarDeclination{NF}()

    # WORK ARRAYS
    initial_time::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
    cos_zenith::Grid = zeros(Grid,nlat_half)
end

function cos_zenith!(
    S::SolarZenithSeason{NF},
    time::DateTime,
    geometry::Geometry
) where NF

    (;sinlat, coslat, lat) = geometry
    (;cos_zenith, length_of_day, length_of_year) = S

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? time : S.initial_time[]
    g = year_angle(NF, time_of_year, length_of_day, length_of_year)
    
    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ,cosδ = sincos(δ)

    local h₀::NF                # hour angle sunrise to sunset
    local cos_zenith_j::NF      # at latitude j

    rings = eachring(cos_zenith)
    @inbounds for (j,ring) in enumerate(rings)
            
        ϕ = lat[j]
        h₀ = abs(δ) + abs(ϕ) < π/2 ?        # polar day/night?
        acos(-tan(ϕ) * tan(δ)) :            # if not: calculate length of day
        ϕ*δ > 0 ? π : 0                     # polar day if signs are equal, otherwise polar night
        
        sinϕ, cosϕ = sinlat[j], coslat[j]
        cos_zenith_j = max(0,h₀*sinδ*sinϕ + cosδ*cosϕ*sin(h₀))
        cos_zenith_j /= π

        for ij in ring
            cos_zenith[ij] = cos_zenith_j
        end
    end
end