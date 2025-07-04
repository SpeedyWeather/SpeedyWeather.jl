# function barrier for GPU
function SpeedyWeather.vorticity_flux_curldiv!(
    diagn::DiagnosticVariables,
    coriolis::SpeedyWeather.AbstractCoriolis,
    geometry::SpeedyWeather.Geometry,
    S::SpeedyWeather.SpectralTransform{NF, <:GPU};
    div::Bool=true,     # also calculate div of vor flux?
    add::Bool=false) where {NF}   # accumulate in vor/div tendencies?
    
    return SpeedyWeather.vorticity_flux_curldiv_kernel!(diagn, coriolis, geometry, S; div, add)
end