"""
Interpolate gridded variable dims: (lonlat,σ) 
to pressure levels: (lonlat,p) with extrapolation
"""
function interpolate_pressure_levels(
    in_field::SpeedyWeather.AbstractField, # input with dims (lonlat, σ-levels)
    p::AbstractVector,
    vars::Variables,
    model::AbstractModel,
)
    @assert issorted(p)

    n_plev = length(p)
    σ      = model.geometry.σ_levels_full
    pₛ     = vars.parameterizations.surface_pressure

    out_field = vars.scratch.grid.a[:,1:n_plev]

    for ij in eachindex(in_field[:,1])
        column_interpolate_pressure_levels!(ij,out_field,in_field,pₛ,σ,p,n_plev)
    end
    
    return out_field
end

"""
Column wise vertical interpolation
"""
function column_interpolate_pressure_levels!(
    ij::Int,
    out_field::SpeedyWeather.AbstractField{T,2}, # input with dims (lonlat, σ-levels)
    in_field::SpeedyWeather.AbstractField{T,2}, # input with dims (lonlat, σ-levels)
    pₛ::SpeedyWeather.AbstractField{T,1},
    σ::AbstractVector{T},
    p::AbstractVector,
    n_plev::Int,
) where T

    for pk in 1:n_plev

        σij = p[pk] / pₛ[ij]

        if σij < σ[1]
            out_field[ij,pk] = in_field[ij,1]

        elseif σij >= σ[end]
            out_field[ij,pk] = in_field[ij,end]

        else
            k = 1
            while (σ[k] < σij) k += 1 end
                    
            Δσ, w1, w2 = interpolation_weights(σij,σ[k-1],σ[k])
            out_field[ij,pk] = ( w1*in_field[ij,k-1] + w2*in_field[ij,k] ) / Δσ

            end
    end
    
    return out_field
end

"""
Calculate weights assuming σ₁ <= σ <= σ₂

returns:

Δσ  = σ₂ - σ₁
Δσ₁ = σ  - σ₁  
Δσ₂ = σ₂ - σ
"""
@inline function interpolation_weights(σ,σ₁,σ₂)
    return σ₂ - σ₁, σ - σ₁, σ₂ - σ
end

"""
Log interpolate gridded variable dims: (lonlat,σ) 
to pressure levels: (lonlat,p) with extrapolation
"""
function log_interpolate_pressure_levels(
    in_field::SpeedyWeather.AbstractField, # input with dims (lonlat, σ-levels)
    p::AbstractVector,
    vars::Variables,
    model::AbstractModel,
)

    @assert issorted(p)
    
    n_plev = length(p)
    σ      = model.geometry.σ_levels_full
    pₛ     = vars.parameterizations.surface_pressure

    out_field = vars.scratch.grid.a[:,1:n_plev]

    for ij in eachindex(in_field[:,1])
        column_log_interpolate_pressure_levels!(ij,out_field,in_field,pₛ,σ,p,n_plev)
    end
    
    return out_field
end

"""
Column wise logarithmic vertical interpolation
"""
function column_log_interpolate_pressure_levels!(
    ij::Int,
    out_field::SpeedyWeather.AbstractField{T,2}, # input with dims (lonlat, σ-levels)
    in_field::SpeedyWeather.AbstractField{T,2}, # input with dims (lonlat, σ-levels)
    pₛ::SpeedyWeather.AbstractField{T,1},
    σ::AbstractVector{T},
    p::AbstractVector,
    n_plev::Int,
) where T

    for pk in 1:n_plev

        σij = p[pk] / pₛ[ij]

        if σij < σ[1]
            out_field[ij,pk] = in_field[ij,1]

        elseif σij >= σ[end]
            out_field[ij,pk] = in_field[ij,end]

        else
            k = 1
            while (σ[k] < σij) k += 1 end
                    
            Δσ, w1, w2 = interpolation_log_weights(σij,σ[k-1],σ[k])
            out_field[ij,pk] = ( w1*in_field[ij,k-1] + w2*in_field[ij,k] ) / Δσ

            end
    end
    
    return out_field
end

"""
Calculate log weights assuming σ₁ <= σ <= σ₂

returns:

Δσ  = log(σ₂/σ₁) = log(σ₂) - log(σ₁)
Δσ₁ = log(σ/σ₁)  = log(σ)  - log(σ₁)
Δσ₂ = log(σ₂/σ)  = log(σ₂) - log(σ)
"""
@inline function interpolation_log_weights(σ,σ₁,σ₂)
    return log(σ₂/σ₁), log(σ/σ₁), log(σ₂/σ)
end