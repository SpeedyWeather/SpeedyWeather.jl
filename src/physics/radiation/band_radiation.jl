
struct BandRadiation{A, C, K, E, M} 
    gases  :: A
    gas_concentration :: C
    gas_absorption_SW :: K
    gas_absorption_LW :: K
    emissivity_SW :: E
    emissivity_LW :: E
    bandmodel :: M
end

function BandRadiation(grid; gases = (:H2O, :CO2, :O3), 
                       bandmodel = GeneralBandModel(limits_LW = ClimLabBandLW(), limits_SW = ClimLabBandSW()))
    
    nlayers   = size(grid)
    nbands_LW = bandsizeLW(bandmodel)
    nbands_SW = bandsizeSW(bandmodel)

    gas_concentration = zeros(eltype(grid), length(gases), nlayers)
    gas_absorption_LW = zeros(eltype(grid), length(gases), nbands_LW)
    gas_absorption_SW = zeros(eltype(grid), length(gases), nbands_SW)

    emissivity_LW     = ones(eltype(grid), nbands_LW, nlayers)
    emissivity_SW     = ones(eltype(grid), nbands_SW, nlayers)
    
    for (idx, gas) in enumerate(gases)
        gas_concentration[idx, :] .= gas_distribution(grid.p_center, gas)
        for band in 1:nbands_LW
            gas_absorption_LW[idx, band] = getabsorption_LW(gas, band, bandmodel)
        end
        for band in 1:nbands_SW
            gas_absorption_SW[idx, band] = getabsorption_SW(gas, band, bandmodel)
        end
    end

    return BandRadiation(gases, gas_concentration, gas_absorption_SW, gas_absorption_LW, emissivity_SW, emissivity_LW, bandmodel)
end

calculate_absorbers!(model, ::GreyRadiation)     = nothing

function calculate_absorbers!(model, radiation)
    for (idx, gas) in enumerate(radiation.gases)
        radiation.gas_concentration[idx, :] .= gas_distribution(model.grid.p_center, gas, T = model.temperature)
    end
end

calculate_emissivities!(model, ::FixedEmissivity) = nothing

function calculate_emissivities!(model, radiation::GreyRadiation) 
    grid = model.grid
    radiation.emissivity   .= 1.0 .- exp.( - grid.Δp_center .* radiation.κ / g) 
    radiation.emissivity[1] = 1.0
    return nothing
end

function calculate_emissivities!(model, radiation::BandRadiation) 

    grid = model.grid

    N = size(grid)

    T   = model.temperature
    nLW = bandsizeLW(radiation.bandmodel)
    nSW = bandsizeSW(radiation.bandmodel)

    for band in 1:nLW
        for i in 2:N
            κ = 0.0
            for gas in 1:length(radiation.gases)
                κ += radiation.gas_concentration[gas, i] * radiation.gas_absorption_LW[gas, band]
            end
            ρ = grid.p_center[i] / (R₀ / M * T[i])
            radiation.emissivity_LW[band, i] = 1 - exp(- κ * ΔpC(grid, i) / g) # ΔhC(grid, i) * ρ) #
        end
    end

    for band in 1:nSW
        for i in 2:N
            κ = 0.0
            for gas in 1:length(radiation.gases)
                κ += radiation.gas_concentration[gas, i] * radiation.gas_absorption_SW[gas, band]
            end
            ρ = grid.p_center[i] / (R₀ / M * T[i])
            radiation.emissivity_SW[band, i] = 1 - exp(- κ * ΔpC(grid, i) / g) #* ΔhC(grid, i) * ρ)
        end
    end

    return nothing
end
