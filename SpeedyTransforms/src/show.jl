# from BenchmarkTools.jl but using KB instead of KiB etc
function prettymemory(b)
    if b < 1.0e3
        return string(b, " bytes")
    elseif b < 1.0e6
        value, units = b * 1.0e-3, "KB"
    elseif b < 1.0e9
        value, units = b * 1.0e-6, "MB"
    else
        value, units = b * 1.0e-9, "GB"
    end
    return string(@sprintf("%.2f", value), " ", units)
end

function Base.show(io::IO, S::SpectralTransform{NF, AR, AT}) where {NF, AR, AT}
    (; spectrum, grid, nlayers, architecture) = S
    (; lmax, mmax) = spectrum   # 1-based max degree/order of harmonics
    (; nlat_half) = grid
    Grid = nonparametric_type(grid)

    # add information about size of Legendre polynomials and scratch memory
    polysize_str = prettymemory(Base.summarysize(S.legendre_polynomials))
    memorysize_str = prettymemory(
        sizeof(S.scratch_memory.north) +      # add all scratch_memories
            sizeof(S.scratch_memory.south) +
            sizeof(S.scratch_memory.column.north) +
            sizeof(S.scratch_memory.column.south)
    )

    dealias = get_dealiasing(mmax - 1, nlat_half) # -1 for zero-based
    truncations = ["<linear", "linear", "quadratic", "cubic", ">cubic"]
    truncation = truncations[clamp(floor(Int, dealias) + 1, 1, 5)]
    dealiasing = @sprintf("%.3g", dealias)

    param_str = "{$NF, $AT}"
    println(io, styled"{warning:SpectralTransform}{note:$param_str}")
    println(io, styled"├ {info:Spectral}:     T$(mmax - 1), $(lmax)x$(mmax) LowerTriangularMatrix{note:\{Complex\{$NF\}\}}")
    println(io, styled"├ {info:Grid}:         Field{note:\{$NF\}}, $(RingGrids.get_nlat(grid))-ring $Grid")
    println(io, styled"├ {info:Truncation}:   dealiasing = $dealiasing {note:($truncation)}")
    println(io, styled"├ {info:Legendre}:     Polynomials $polysize_str, shortcut: $(short_name(S.LegendreShortcut))")
    println(io, styled"├ {info:Architecture}: $architecture")
    print(io, styled"└ {info:Memory}:       for $nlayers layers {note:($memorysize_str)}")
    return nothing
end

function Base.summary(io::IO, M::ScratchMemory)
    AT = nonparametric_type(typeof(M.north))
    NF = eltype(M.north)
    return print(io, Base.dims2string((size(M.north)..., 2)), " ScratchMemory{$AT{$NF,...}, ...}")
end

function Base.show(io::IO, S::MatrixSpectralTransform{NF, AR, AT}) where {NF, AR, AT}
    (; spectrum, grid, nlayers, architecture) = S
    (; lmax, mmax) = spectrum   # 1-based max degree/order of harmonics
    (; nlat_half) = grid
    Grid = nonparametric_type(grid)

    # add information about size of Legendre polynomials and scratch memory
    matrixsize_str = prettymemory(2 * Base.sizeof(S.forward))
    scratchsize_str = prettymemory(Base.sizeof(S.scratch_memory))

    dealias = get_dealiasing(mmax - 1, nlat_half) # -1 for zero-based
    truncations = ["<linear", "linear", "quadratic", "cubic", ">cubic"]
    truncation = truncations[clamp(floor(Int, dealias) + 1, 1, 5)]
    dealiasing = @sprintf("%.3g", dealias)

    param_str = "{$NF, $AT}"
    println(io, styled"{warning:MatrixSpectralTransform}{note:$param_str}")
    println(io, styled"├ {info:Spectral}:     T$(mmax - 1), $(lmax)x$(mmax) LowerTriangularMatrix{note:\{Complex\{$NF\}\}}")
    println(io, styled"├ {info:Grid}:         Field{note:\{$NF\}}, $(RingGrids.get_nlat(grid))-ring $Grid")
    println(io, styled"├ {info:Truncation}:   dealiasing = $dealiasing {note:($truncation)}")
    println(io, styled"├ {info:Architecture}: $architecture")
    print(io, styled"└ {info:Memory}:       $matrixsize_str matrices, $scratchsize_str scratch ($nlayers layers)")
    return nothing
end
