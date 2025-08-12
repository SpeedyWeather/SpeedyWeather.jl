"""Legendre shortcut is the truncation of the loop over the order m of the spherical harmonics
in the Legendre transform. For reduced grids with as few as 4 longitudes around the poles
(HEALPix grids) or 20 (octahedral grids) the higher wavenumbers in large orders m do not
project (significantly) onto such few longitudes. For performance reasons the loop over m can
therefore but truncated early. A Legendre shortcut `<: AbstractLegendreShortcut`
is implemented as a functor that returns the 0-based maximum order m to retain per latitude ring,
i.e. to be used for m in 0:mmax_truncation.

New shortcuts can be added by defining `struct LegendreShortcutNew <: AbstractLegendreShortcut end`
and the functor method `LegendreShortcutNew(nlon::Integer, lat::Real) = ...`, with `nlon` the
number of longitude points on that ring, and `latd` its latitude in degrees (-90˚ to 90˚N).
Many implementations may not use the latitude `latd` but it is included for compatibility.
If unused set to default value to 0. Also define `short_name(::Type{<:LegendreShortcutNew}) = "new"`.

Implementions are `LegendreShortcutLinear`, `LegendreShortcutQuadratic`, `LegendreShortcutCubic`,
`LegendreShortcutLinQuadCosLat²` and `LegendreShortcutLinCubCoslat`."""
abstract type AbstractLegendreShortcut end
short_name(s::Type{<:AbstractLegendreShortcut}) = string(s)
short_name(s::AbstractLegendreShortcut) = short_name(typeof(s))

struct LegendreShortcutLinear <: AbstractLegendreShortcut end
"""$(TYPEDSIGNATURES)
Linear Legendre shortcut, truncates the Legendre loop over order m to nlon/2."""
LegendreShortcutLinear(nlon::Integer, latd::Real=0) = nlon÷2
short_name(::Type{<:LegendreShortcutLinear}) = "linear"

struct LegendreShortcutQuadratic <: AbstractLegendreShortcut end
"""$(TYPEDSIGNATURES)
Quadratic Legendre shortcut, truncates the Legendre loop over order m to nlon/3."""
LegendreShortcutQuadratic(nlon::Integer, latd::Real=0) = nlon÷3
short_name(::Type{<:LegendreShortcutQuadratic}) = "quadratic"

struct LegendreShortcutCubic <: AbstractLegendreShortcut end
"""$(TYPEDSIGNATURES)
Cubic Legendre shortcut, truncates the Legendre loop over order m to nlon/4."""
LegendreShortcutCubic(nlon::Integer, latd::Real=0) = nlon÷4
short_name(::Type{<:LegendreShortcutCubic}) = "cubic"

struct LegendreShortcutLinQuadCoslat² <: AbstractLegendreShortcut end
"""$(TYPEDSIGNATURES)
Linear-Quadratic Legendre shortcut, truncates the Legendre loop over order m to nlon/(2 + cosd(latd)^2)."""
LegendreShortcutLinQuadCoslat²(nlon::Integer, latd::Real) = floor(Int, nlon/(2 + cosd(latd)^2))
short_name(::Type{<:LegendreShortcutLinQuadCoslat²}) = "linquadcoslat²"

struct LegendreShortcutLinCubCoslat <: AbstractLegendreShortcut end
"""$(TYPEDSIGNATURES)
Linear-Cubic Legendre shortcut, truncates the Legendre loop over order m to nlon/(2 + 2cosd(latd))."""
LegendreShortcutLinCubCoslat(nlon::Integer, latd::Real) = floor(Int, nlon/(2 + 2cosd(latd)))
short_name(::Type{<:LegendreShortcutLinCubCoslat}) = "lincubcoslat"