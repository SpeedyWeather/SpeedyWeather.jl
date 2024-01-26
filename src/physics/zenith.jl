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

"""Struct containing all parameters for the calculation of the
cos of the solar zenith angle, cos_zenith."""
Base.@kwdef struct SolarZenithAngle{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractZenith{NF,Grid}
    # DIMENSIONS
    nlat_half::Int

    # OPTIONS
    daily_cycle::Bool = true
    length_of_day::Second = Hour(24)
    seasonal_cycle::Bool = true
    length_of_year::Second = Day(365.25)

    # COEFFICIENTS
    solar_declination::SinSolarDeclination{NF} = SinSolarDeclination{NF}()
    time_correction::SolarTimeCorrection{NF} = SolarTimeCorrection{NF}()

    # WORK ARRAYS
    initial_time::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
    cos_zenith::Grid = zeros(Grid,nlat_half)
end

# generator function
function SolarZenithAngle(SG::SpectralGrid,P::AbstractPlanet;kwargs...)
    (;NF, Grid, nlat_half) = SG
    (;daily_cycle, seasonal_cycle, length_of_day, length_of_year) = P
    solar_declination = SinSolarDeclination(SG,P)
    SolarZenithAngle{NF,Grid{NF}}(;nlat_half,daily_cycle,seasonal_cycle,
        length_of_day,length_of_year,solar_declination,kwargs...)
end

function Base.show(io::IO,L::AbstractZenith)
    println(io,"$(typeof(L)) <: AbstractZenith")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

function initialize!(
    S::SolarZenithAngle,
    initial_time::DateTime,
    model::ModelSetup
)
    S.initial_time[] = initial_time
    cos_zenith!(S,initial_time,model.geometry)
end

function cos_zenith!(
    S::SolarZenithAngle{NF},
    time::DateTime,
    geometry::Geometry
) where NF

    (;sinlat, coslat, lat, lons) = geometry
    (;cos_zenith, length_of_day, length_of_year) = S

    # convert day and year in seconds to radians
    day2rad = convert(NF,2π/length_of_day.value)
    year2rad = convert(NF,2π/length_of_year.value)

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    # but use the initial time in case of no seasonal cycle
    time_of_year = S.seasonal_cycle ? time : S.initial_time[]
    sec_of_day = Second(Dates.Time(time_of_year).instant).value
    g = year2rad*(Dates.dayofyear(time_of_year)*length_of_day.value + sec_of_day)

    # time correction [radians] due to the equation of time (sunrise/set oscillation)
    tc = S.time_correction(g)
    
    # solar hour angle at 0˚E (longtiude offset added later)
    noon_in_sec = length_of_day.value ÷ 2
    sec_of_day = Second(Dates.Time(time).instant).value
    solar_hour_angle_0E = (sec_of_day - noon_in_sec)*day2rad + tc

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ,cosδ = sincos(δ)

    rings = eachring(cos_zenith)

    if S.daily_cycle
        for (j,ring) in enumerate(rings)                         
            sinϕ, cosϕ = sinlat[j], coslat[j]            # sin, cos of latitude
            for ij in ring
                h = solar_hour_angle_0E + lons[ij]      # solar hour angle at longitude λ in radians
                cos_zenith[ij] = max(0,sinδ*sinϕ + cosδ*cosϕ*cos(h))
            end
        end
    else
        for (j,ring) in enumerate(rings)
            
            ϕ = lat[j]
            h₀ = abs(δ) + abs(ϕ) < π/2 ?    # there is a sunset / sunrise
            acos(-tan(ϕ) * tan(δ)) :
            ϕ * δ > 0 ? π : 0.0             # polar day or polar night
            sinh₀_h₀ = h₀ == 0 ? 1 : sin(h₀)/h₀
            
            sinϕ, cosϕ = sinlat[j], coslat[j]            # sin, cos of latitude
            cos_zenith_j = max(0,sinδ*sinϕ + cosδ*cosϕ*sinh₀_h₀)

            for ij in ring
                cos_zenith[ij] = cos_zenith_j
            end
        end
    end
end