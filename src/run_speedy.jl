function run_speedy(::Type{T}=Float64;      # number format
                    kwargs...               # all additional parameters
                    ) where {T<:AbstractFloat}

    P = Params(T=T,kwargs...)
    C = Constants{T}(P)
    G = GeoSpectral{T}(P)
    B = Boundaries{T}(P,G)

    Prog = initial_conditions(P,B,G)

   

    # # TODO
    # Prog = PrognosticVars{T}()
    # Diag = DiagnosticVars{T}()

    #time_stepping!()

    return G,B
end


println("hello world")