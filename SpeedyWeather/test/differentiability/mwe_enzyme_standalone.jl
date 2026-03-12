## Standalone Enzyme MWE: GC-pointer phi nodes in copy! via Dict iteration
##
## Background: SpeedyWeather's Base.copy!(::PrognosticVariables, ...) iterates
## over Dict-like fields using `for (key, value) in pairs(...)`. On GitHub
## Actions (ubuntu-latest, AMD EPYC, Julia 1.10) this occasionally causes
## an Enzyme internal error ("Could not analyze GC behavior" / nodecayed_phis!).
## This MWE attempts to reproduce the pattern without SpeedyWeather.
##
## Run: julia +1.10 mwe_enzyme_standalone.jl
## Expected on GH Actions: may fail with EnzymeInternalError

using Enzyme

Enzyme.Compiler.VERBOSE_ERRORS[] = true

# ============================================================
# Minimal struct mimicking SpeedyWeather's prognostic variables
# with Dict-like ocean/land fields
# ============================================================

mutable struct FieldContainer
    data::Dict{Symbol, Vector{Float64}}
end

mutable struct StateVec
    ocean::FieldContainer
    land::FieldContainer
    core::Vector{Float64}
end

function Base.copy!(dst::StateVec, src::StateVec)
    # Pattern 1: loop over Dict entries (suspected culprit)
    for (key, value) in pairs(src.ocean.data)
        dst.ocean.data[key] .= value
    end
    for (key, value) in pairs(src.land.data)
        dst.land.data[key] .= value
    end
    # Pattern 2: direct vector copy (baseline, expected to work)
    dst.core .= src.core
    return dst
end

function f_copy!(dst, src)
    copy!(dst, src)
    return nothing
end

n = 8
src = StateVec(
    FieldContainer(Dict(:sst => rand(n), :ice => rand(n), :flux => rand(n))),
    FieldContainer(Dict(:soil => rand(n), :temp => rand(n))),
    rand(n),
)
dst = StateVec(
    FieldContainer(Dict(:sst => zeros(n), :ice => zeros(n), :flux => zeros(n))),
    FieldContainer(Dict(:soil => zeros(n), :temp => zeros(n))),
    zeros(n),
)
ddst = StateVec(
    FieldContainer(Dict(:sst => ones(n), :ice => ones(n), :flux => ones(n))),
    FieldContainer(Dict(:soil => ones(n), :temp => ones(n))),
    ones(n),
)
dsrc = StateVec(
    FieldContainer(Dict(:sst => zeros(n), :ice => zeros(n), :flux => zeros(n))),
    FieldContainer(Dict(:soil => zeros(n), :temp => zeros(n))),
    zeros(n),
)

println("Testing Enzyme autodiff through copy! with Dict iteration...")
flush(stdout)

try
    autodiff(Reverse, f_copy!, Const,
        Duplicated(dst, ddst),
        Duplicated(src, dsrc),
    )
    println("SUCCESS: autodiff through Dict-iterating copy! worked")
catch e
    println("FAILED with $(typeof(e)):")
    println(e)
end

# ============================================================
# Variant: only the Dict loop part (no core vector)
# ============================================================

mutable struct DictState
    data::Dict{Symbol, Vector{Float64}}
end

function g_copy!(dst::DictState, src::DictState)
    for (key, value) in pairs(src.data)
        dst.data[key] .= value
    end
    return nothing
end

src2 = DictState(Dict(:a => rand(n), :b => rand(n), :c => rand(n)))
dst2 = DictState(Dict(:a => zeros(n), :b => zeros(n), :c => zeros(n)))
ddst2 = DictState(Dict(:a => ones(n), :b => ones(n), :c => ones(n)))
dsrc2 = DictState(Dict(:a => zeros(n), :b => zeros(n), :c => zeros(n)))

println()
println("Testing Enzyme autodiff through bare Dict-loop copy!...")
flush(stdout)

try
    autodiff(Reverse, g_copy!, Const,
        Duplicated(dst2, ddst2),
        Duplicated(src2, dsrc2),
    )
    println("SUCCESS: bare Dict-loop copy! worked")
catch e
    println("FAILED with $(typeof(e)):")
    println(e)
end
