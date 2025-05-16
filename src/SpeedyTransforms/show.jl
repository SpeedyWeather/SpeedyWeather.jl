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
    (; spectrum, Grid, nlat_half, nlayers) = S
    (; lmax, mmax) = spectrum

    # add information about size of Legendre polynomials and scratch memory
    polysize_str = prettymemory(Base.summarysize(S.legendre_polynomials))
    memorysize_str = prettymemory(
                Base.summarysize(S.scratch_memory_north) +      # add all scratch_memories
                Base.summarysize(S.scratch_memory_south) + 
                Base.summarysize(S.scratch_memory_grid) + 
                Base.summarysize(S.scratch_memory_spec)
            )

    dealias = get_dealiasing(mmax, nlat_half)
    truncations = ["<linear", "linear", "quadratic", "cubic", ">cubic"]
    truncation = truncations[clamp(floor(Int, dealias)+1, 1, 5)]
    dealiasing = @sprintf("%.3g", dealias)

    println(io, "SpectralTransform{$NF, $ArrayType}:")
    println(io, "├ Spectral:   T$(mmax-1), $(lmax)x$(mmax) LowerTriangularMatrix{Complex{$NF}}")
    println(io, "├ Grid:       $(RingGrids.get_nlat(Grid, nlat_half))-ring $Grid{$NF}")
    println(io, "├ Truncation: dealiasing = $dealiasing ($truncation)")
    println(io, "├ Legendre:   Polynomials $polysize_str, shortcut: $(short_name(S.LegendreShortcut))")
    println(io, "└ Memory:     for $nlayers layers ($memorysize_str)")
end