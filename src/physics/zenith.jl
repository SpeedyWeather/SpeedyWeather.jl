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

Base.@kwdef struct SolarTimeCorrection{NF} <: AbstractSolarTimeCorrection{NF}
    a::NF =   0.004297      # the offset +
    s1::NF = -1.837877      # s1*sin(g) +
    c1::NF =  0.107029      # c1*cos(g) +
    s2::NF = -2.340475      # s2*sin(2g) +
    c2::NF = -0.837378      # c2*cos(2g)
end

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

Base.@kwdef struct SolarZenithAngle{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractZenith{NF,Grid}
    # DIMENSIONS
    nlat_half::Int

    # OPTIONS
    daily_cycle::Bool = true
    length_of_day::Second = Hour(24)
    seasonal_cycle::Bool = true
    length_of_year::Second = Day(365.25)

    # COEFFICIENTS
    frac_year2rad::NF = 2π/length_of_year.value
    solar_declination::SolarDeclination{NF} = SolarDeclination{NF}()
    time_correction::SolarTimeCorrection{NF} = SolarTimeCorrection{NF}()

    # WORK ARRAYS
    cos_zenith::Grid = zeros(Grid,nlat_half)
end

# generator function
function SolarZenithAngle(SG::SpectralGrid,P::AbstractPlanet;kwargs...)
    (;NF, Grid, nlat_half) = SG
    (;daily_cycle, seasonal_cycle) = P
    SolarZenithAngle{NF,Grid{NF}}(;nlat_half,daily_cycle,seasonal_cycle,kwargs...)
end

function Base.show(io::IO,L::AbstractZenith)
    println(io,"$(typeof(L)) <: AbstractZenith")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

function cos_zenith!(
    S::SolarZenithAngle,
    time::DateTime,
    geometry::Geometry
)

    (;sinlat, coslat, lons) = geometry
    (;cos_zenith, length_of_day) = S

    # angular fraction of year [0...2π] for Jan-01 to Dec-31
    sec_of_day = Second(Dates.Time(time).instant).value
    g = S.frac_year2rad*(Dates.dayofyear(time)*length_of_day.value + sec_of_day)
    tc = S.time_correction(g)
    noon_in_sec = length_of_day.value ÷ 2
    solar_hour_angle_0E = (sec_of_day - noon_in_sec)*S.frac_year2rad + tc
    δ = S.solar_declination(g)
    sinδ,cosδ = sincos(δ)

    rings = eachring(cos_zenith)
    for (j,ring) in enumerate(rings)                         
        sinϕ, cosϕ = sinlat[j], coslat[j]            # sin, cos of latitude
        for ij in ring
            h = solar_hour_angle_0E + lons[ij]      # solar hour angle at longitude λ in radians
            cos_zenith[ij] = max(0,sinδ*sinϕ + cosδ*cosϕ*cos(h))
        end
    end
end