# stolen from BenchmarkTools.jl
function prettymemory(b)
  if b < 1024
      return string(b, " bytes")
  elseif b < 1024^2
      value, units = b / 1024, "KiB"
  elseif b < 1024^3
      value, units = b / 1024^2, "MiB"
  else
      value, units = b / 1024^3, "GiB"
  end
  return string(@sprintf("%.2f", value), " ", units)
end

function Base.show(io::IO,S::SpectralTransform{NF}) where NF
    (;lmax,mmax,Grid,nlat_half) = S

    # add information about size of Legendre polynomials
    s = S.recompute_legendre ? Base.summarysize(S.Λ) : Base.summarysize(S.Λs)
    s_str = prettymemory(s)

    println(io,"$(typeof(S))(")
    println(io,"  Spectral: T$mmax, $(lmax+1)x$(mmax+1) LowerTriangularMatrix{Complex{$NF}}")
    println(io,"  Grid:     $(RingGrids.get_nlat(Grid,nlat_half))-ring $Grid{$NF}")
      print(io,"  Legendre: recompute polynomials $(S.recompute_legendre), $s_str)")
end
