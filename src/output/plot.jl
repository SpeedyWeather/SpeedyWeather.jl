# to avoid ambiguity with exporting plot
plot(L::LowerTriangularMatrix{T}; kwargs...) where T = LowerTriangularMatrices.plot(L;kwargs...)
plot(G::AbstractGrid{T}; kwargs...) where T = RingGrids.plot(G;kwargs...)