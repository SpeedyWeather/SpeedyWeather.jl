struct NoForcing{NF} <: AbstractForcing{NF} end
NoForcing(SG::SpectralGrid) = NoForcing{SG.NF}()

function initialize!(   forcing::NoForcing,
                        model::ModelSetup)
    return nothing
end