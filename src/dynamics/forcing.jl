struct NoForcing <: AbstractForcing end

function initialize!(   forcing::AbstractForcing,
                        model::ModelSetup)
    return nothing
end