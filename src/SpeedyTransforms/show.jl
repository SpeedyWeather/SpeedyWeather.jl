# from BenchmarkTools.jl but using KB instead of KiB etc
function prettymemory(b)
    if b < 1e3
        return string(b, " bytes")
    elseif b < 1e6
        value, units = b * 1e-3, "KB"
    elseif b < 1e9
        value, units = b * 1e-6, "MB"
    else
        value, units = b * 1e-9, "GB"
    end
    return string(@sprintf("%.2f", value), " ", units)
end

function Base.show(io::IO, S::SpectralTransform{NF, ArrayType}) where {NF, ArrayType}
    (; lmax, mmax, Grid, nlat_half, nlayers) = S

    # add information about size of Legendre polynomials
    s = S.recompute_legendre ? Base.summarysize(S.Λ) : Base.summarysize(S.Λs)
    s_str = prettymemory(s)
    m_str = prettymemory(Base.summarysize(S.scratch_memory_north)*2)

    dealias = get_dealiasing(mmax, nlat_half)
    truncations = ["<linear", "linear", "quadratic", "cubic", ">cubic"]
    truncation = truncations[clamp(floor(Int, dealias)+1, 1, 5)]
    dealiasing = @sprintf("%.3g", dealias)

    println(io, "SpectralTransform{$NF, $ArrayType}:")
    println(io, "├ Spectral:   T$mmax, $(lmax+1)x$(mmax+1) LowerTriangularMatrix{Complex{$NF}}")
    println(io, "├ Grid:       $(RingGrids.get_nlat(Grid, nlat_half))-ring $Grid{$NF}")
    println(io, "├ Truncation: dealiasing = $dealiasing ($truncation)")
    println(io, "├ Legendre:   recompute polynomials $(S.recompute_legendre) ($s_str)")
      print(io, "└ Memory:     scratch for $nlayers layers ($m_str)")
end
