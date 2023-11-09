
struct NoRadiation{NF} <: AbstractRadiationModel{NF}

struct Radiation{NF, E, A, S, V} <: AbstractRadiationModel{NF}
    emissivity             :: E
    albedo                 :: A
    solar_forcing          :: S
    upward_flux            :: V
    downward_flux          :: V
    absorbed_upward_flux   :: V
    absorbed_downward_flux :: V
    skip_time_steps        :: Int
end

"""
    NoRadiation(args...) 
    
fallback for models without radiation
"""
NoRadiation(args...) = NoRadiation()

initialize!(radiation::NoRadiation) = radiation

radiation!(column,model,::NoRadiation) = nothing

function Radiation(spectral_grid::SpectralGrid;  
                   solar_forcing = AnnualAverageSolarForcing(spectral_grid), 
                   albedo     = 0.2985, 
                   emissivity = GrayEmissivity(spectral_grid; ϵ = 0.8))

    upward_flux   = similar(temperatures)
    downward_flux = similar(temperatures)
    absorbed_upward_flux   = similar(temperatures)
    absorbed_downward_flux = similar(temperatures)
   
    types = (spectral_grid.NF, typeof(emissivity), typeof(albedo), typeof(solar_forcing), typeof(upward_flux), Int)

    return Radiation{types...}(emissivity,
                               albedo,
                               solar_forcing,
                               upward_flux, 
                               downward_flux, 
                               absorbed_upward_flux, 
                               absorbed_downward_flux, 
                               0)
end

struct GrayEmissivity{NF, E} 
    ϵ :: E
end

"""
    GrayEmissivity(spectral_grid; ϵ = 1.0)

A radiation model with a gray absorption coefficient for infrared equal to
κ m⁻¹ and no absorption in the UV and visible part of the spectrum
"""
GrayEmissivity(spectral_grid::SpectralGrid; ϵ = 1.0) = GreyRadiation(convert(spectral_grid.NF, ϵ))

const GrayRadiation = Radiation{<:Any, <:GrayEmissivity}

emitted_radiation(r::GrayRadiation, T, k) = 2 * emissivity(r, T) * σ * T[k]^4

function upward_flux!(model, radiation::GreyRadiation)
    
    ε   = radiation.emissivity
    T   = model.temperature
    UF  = model.upward_flux
    AUF = model.absorbed_upward_flux

    fill!(UF, 0.0)
    
    N = size(model)
    @inbounds for i in 2:N
        transmitted_radiation = (1 - ε[i-1]) * UF[i-1]
        direct_radiation      = ε[i-1] * σ * T[i-1]^4

        UF[i]  = ( direct_radiation + transmitted_radiation )
        AUF[i] = ε[i] * UF[i] 
    end
end

function downward_flux!(model, radiation::GlobalRadiation)

    ε   = radiation.emissivity
    T   = model.temperature
    DF  = model.downward_flux
    ADF = model.absorbed_downward_flux

    fill!(DF, 0.0)

    N = size(model)
    @inbounds for i in N-1:-1:1
        transmitted_radiation = (1 - ε[i+1]) * DF[i+1]
        direct_radiation      = ε[i+1] * σ * T[i+1]^4

        DF[i]  = ( direct_radiation + transmitted_radiation )
        ADF[i] = ε[i] * DF[i] 
    end

    DF[1]  += (1 - model.albedo) * model.forcing
    ADF[1] += (1 - model.albedo) * model.forcing
end

function emitted_radiation(r::BandRadiation, T, i) 
    E = 0
    for band in 1:bandsizeLW(r.bandmodel)
        E += 2 * r.emissivity_LW[band, i] * blackbodyfraction(band, T[i], r.bandmodel) * σ * T[i]^4
    end
    return E
end

function upward_flux!(model, radiation::BandRadiation; calc_OLR = false)
    
    grid = model.grid

    ε   = radiation.emissivity_LW
    b   = radiation.bandmodel

    T   = model.temperature
    UF  = model.upward_flux
    AUF = model.absorbed_upward_flux 
    
    nLW = bandsizeLW(radiation.bandmodel)

    N = size(model)

    fill!(UF,  0.0)
    fill!(AUF, 0.0)

    UFtemp = zeros(eltype(grid), N)

    olr = zeros(nLW)

    @inbounds for band in 1:nLW
        fill!(UFtemp, 0.0)
        for i in 2:N
            transmitted_radiation = (1.0 - ε[band, i-1]) * UFtemp[i-1]
            direct_radiation      = ε[band, i-1] * blackbodyfraction(band, T[i-1], b) * σ * T[i-1]^4
            UFtemp[i] = ( transmitted_radiation + direct_radiation )
        end
        UF  .+= UFtemp
        AUF .+= UFtemp .* ε[band, :]
        if calc_OLR
            olr[band] = (1 - ε[band, N]) * UFtemp[N] + ε[band, N] *  blackbodyfraction(band, T[N], b) * σ * T[end]^4
        end
    end
    if calc_OLR
        return olr
    end
end

# In the downward flux there is also the SW radiation to consider
function downward_flux!(model, radiation::BandRadiation; calc_ASR = false)

    grid = model.grid
    
    εSW  = radiation.emissivity_SW
    εLW  = radiation.emissivity_LW
    b    = radiation.bandmodel

    T   = model.temperature
    DF  = model.downward_flux
    ADF = model.absorbed_downward_flux

    nLW = bandsizeLW(radiation.bandmodel)
    nSW = bandsizeSW(radiation.bandmodel)

    N = size(model)

    fill!(DF,  0.0)
    fill!(ADF, 0.0)

    DFtemp = zeros(eltype(grid), N)
    
    asr = zeros(nSW)

    # Short wave radiation (incoming)
    @inbounds for band in 1:nSW
        fill!(DFtemp, 0.0)
        DFtemp[N] = model.forcing * fluxfraction(band, b) 

        for i in N-1:-1:2
            transmitted_radiation = (1 - εSW[band, i+1]) * DFtemp[i+1]

            DFtemp[i] = (transmitted_radiation)
            if calc_ASR
                asr[band] += εSW[band, i] * transmitted_radiation
            end
        end

        DFtemp[1] = (1 - εSW[band, 2]) * DFtemp[2] * (1.0 - model.albedo)

        if calc_ASR
            asr[band] += DFtemp[1]
        end

        DF  .+= DFtemp
        ADF .+= DFtemp .* εSW[band, :]
    end

    if calc_ASR
        return asr
    end

    # Long wave radiation (emitted)
    @inbounds for band in 1:nLW
        fill!(DFtemp, 0.0)

        for i in N-1:-1:1
            transmitted_radiation = (1 - εLW[band, i+1]) * DFtemp[i+1]
            direct_radiation      = εLW[band, i+1] * blackbodyfraction(band, T[i+1], b) * σ * T[i+1]^4

            DFtemp[i] = ( direct_radiation + transmitted_radiation )
        end
        DF  .+= DFtemp
        ADF .+= DFtemp .* εLW[band, :]
    end
end

OLR(model, radiation::GlobalRadiation) = (1 - radiation.emissivity[end]) * model.upward_flux[end] + radiation.emissivity[end] * σ * model.temperature[end]^4
OLR(model, radiation::BandRadiation)   = upward_flux!(model, radiation; calc_OLR = true)

ASR(model, radiation::GlobalRadiation) = (1 - model.albedo) * model.forcing
ASR(model, radiation::BandRadiation)   = downward_flux!(model, radiation, calc_ASR = true)