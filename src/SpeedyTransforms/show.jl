function Base.show(io::IO,S::SpectralTransform{NF}) where NF
    (;lmax,Grid,nlat_half) = S
    println(io,"$(typeof(S))(")
    println(io,"  Spectral: T$lmax LowerTriangularMatrix{Complex{$NF}}")
    println(io,"  Grid:     $(RingGrids.get_nlat(Grid,nlat_half))-ring $Grid{$NF}")
      print(io,"  Legendre: recompute polynomials $(S.recompute_legendre))")
end
