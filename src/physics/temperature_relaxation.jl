"""NoBoundaryLayer scheme just passes."""
function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::NoTemperatureRelaxation,
                                    model::PrimitiveEquation)
    return nothing
end

"""NoBoundaryLayer scheme does not need any initialisation."""
function initialise_temperature_relaxation!(K::ParameterizationConstants,
                                            scheme::NoTemperatureRelaxation,
                                            P::Parameters,
                                            G::Geometry)
    return nothing
end 

function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::HeldSuarez,
                                    model::PrimitiveEquation)
    (;temp,temp_tend,ϕj) = column
    (;temp_relax_freq,temp_equilibrium) = model.parameterization_constants

    @inbounds for k in eachlayer(column)
        kₜ = temp_relax_freq[k,ϕj]              # (inverse) relaxation time scale
        Teq = temp_equilibrium[k,ϕj]            # eq. temperature [K] at level k, latitude ring ϕj
        temp_tend[k] -= kₜ*(temp[k] - Teq)      # Held and Suarez 1996, equation 2
    end
end
        