# Revert FFTA docs compat restriction

> Status: **completed**. FFTA dependency and compat bound removed from `docs/Project.toml`.

Date of initial draft: 2026-11-06

Base revision: 661e2a34feeb038435abf3870c0f02f54a4a4375

## Originating prompt

> Can you revert PR1372? The upstream issue https://github.com/JuliaMath/FFTA.jl/issues/146 has been merged and FFTA v0.3.3 has been released, so hopefully to version constraint isn't needed anymore, open as PR

## Revision log

## Problem description

#1372 added FFTA as a direct docs dependency only to cap it at 0.3.0 - 0.3.1, because FFTA 0.3.2's
`plan_rfft` was ambiguous with FFTW's. FFTA 0.3.3 fixes this (FFTA.jl#146).

## Summary of changes

Remove FFTA from `[deps]` and `[compat]` in `docs/Project.toml`; it is again only pulled in
transitively (Makie → KernelDensity) and resolves to ≥ 0.3.3.

## Testing and verification

Docs CI build on this PR.
