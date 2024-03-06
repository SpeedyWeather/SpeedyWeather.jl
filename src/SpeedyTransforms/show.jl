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

function Base.show(io::IO, S::SpectralTransform{NF}) where NF
    (;lmax, mmax, Grid, nlat_half) = S

    # add information about size of Legendre polynomials
    s = S.recompute_legendre ? Base.summarysize(S.Λ) : Base.summarysize(S.Λs)
    s_str = prettymemory(s)

    dealias = get_dealiasing(mmax, nlat_half)
    truncations = ["<linear", "linear", "quadratic", "cubic", ">cubic"]
    truncation = truncations[clamp(floor(Int, dealias)+1, 1, 5)]
    dealiasing = @sprintf("%.3g", dealias)

    println(io, "$(typeof(S)):")
    println(io, "├ Spectral:   T$mmax, $(lmax+1)x$(mmax+1) LowerTriangularMatrix{Complex{$NF}}")
    println(io, "├ Grid:       $(RingGrids.get_nlat(Grid, nlat_half))-ring $Grid{$NF}")
    println(io, "├ Truncation: dealiasing = $dealiasing ($truncation)")
      print(io, "└ Legendre:   recompute polynomials $(S.recompute_legendre) ($s_str)")
end
