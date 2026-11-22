abstract type AbstractOceanHeatTransport <: AbstractModelComponent end

export NoOceanHeatTransport, OceanHeatDiffusion

"""No horizontal heat transport in the slab ocean mixed layer (default). The
sea surface temperature tendency of [`SlabOcean`](@ref) then only consists of the surface fluxes."""
struct NoOceanHeatTransport <: AbstractOceanHeatTransport end

initialize!(::NoOceanHeatTransport, ::AbstractOcean, ::AbstractModel) = nothing
heat_transport!(dsst, vars, ::NoOceanHeatTransport, ocean, model) = nothing

"""Horizontal diffusion of the mixed-layer temperature of a [`SlabOcean`](@ref), a simple
parameterization of ocean heat transport by eddies and currents. Conservative finite-volume
discretization on ring grids, precomputed at `initialize!`. Heat (not temperature) is the conserved
quantity, there is no flux through coastlines (every face is weighted by the smaller ocean fraction
of its two cells) and the transport continues under sea ice. Between neighbouring rings the temperature
of both rings is linearly interpolated in longitude to the midpoints of the segments of the
shared face. To be used like

    ocean = SlabOcean(spectral_grid, heat_transport = OceanHeatDiffusion(spectral_grid))

The transport is set by the depth-integrated diffusivity `diffusivity_depth` [m³/s], the flux is
`density * specific_heat_capacity * diffusivity_depth * ∇T`, so the transport does not depend on the mixed-layer
depth `H` and the diffusivity is `D = diffusivity_depth / H`. The default 1.25e6 m³/s is Codron's (2012, Clim. Dyn.,
doi:10.1007/s00382-011-1031-3) 25000 m²/s for a 50 m layer. Alternatively set `diffusivity` [m²/s] directly, independent
of the depth. Both are independent of the planet's radius. The land-sea mask has to be set before `initialize!`.
Any grid with iso-latitude rings and equidistant longitudes per ring is supported, rings at the poles are not.
The scheme is explicit in time, so at `initialize!` the number of substeps per ocean time step is chosen such that
the stability number is below 0.5. Fields and options are $(TYPEDFIELDS)"""
@kwdef mutable struct OceanHeatDiffusion{NF, MatrixInt, MatrixNF, VectorNF} <: AbstractOceanHeatTransport
    "[OPTION] Depth-integrated diffusivity, diffusivity times mixed-layer depth [m³/s]"
    diffusivity_depth::NF = 1.25e6

    "[OPTION] Diffusivity [m²/s], overrides `diffusivity_depth` if not `NaN`"
    diffusivity::NF = NaN

    "[OPTION] Largest stability number (stable below 1) allowed per substep"
    max_stability::NF = 0.5

    "[DERIVED] Diffusivity used [m²/s]"
    diffusivity_used::NF = NaN

    "[DERIVED] Number of explicit substeps per ocean time step"
    substeps::Int = 1

    "[DERIVED] Stencil: linear indices of the grid points each grid point's tendency depends on, padded with itself"
    columns::MatrixInt = zeros(Int, 0, 0)

    "[DERIVED] Stencil: weights for `columns` [1/s], zero for padding, rows sum to zero"
    weights::MatrixNF = zeros(NF, 0, 0)

    "[DERIVED] Area of every grid cell [m²]"
    area::VectorNF = zeros(NF, 0)

    "[DERIVED] Largest sum of absolute stencil weights [1/s], bounds the largest eigenvalue"
    max_rate::NF = 0

    "[DERIVED] Work arrays for the temperature increment of the substeps"
    increment::VectorNF = zeros(NF, 0)
    increment_next::VectorNF = zeros(NF, 0)
end

"""$(TYPEDSIGNATURES) Generator using the number format, array type and architecture of `spectral_grid`."""
function OceanHeatDiffusion(SG::SpectralGrid; kwargs...)
    (; NF, architecture) = SG
    columns = on_architecture(architecture, zeros(Int, 0, 0))
    weights = on_architecture(architecture, zeros(NF, 0, 0))
    area = on_architecture(architecture, zeros(NF, 0))
    return OceanHeatDiffusion{NF, typeof(columns), typeof(weights), typeof(area)}(;
        columns, weights, area, increment = copy(area), increment_next = copy(area), kwargs...
    )
end

function initialize!(heat_transport::OceanHeatDiffusion{NF}, ocean::AbstractOcean, model::AbstractModel) where {NF}
    D = isnan(heat_transport.diffusivity) ? heat_transport.diffusivity_depth / ocean.mixed_layer_depth : heat_transport.diffusivity
    heat_transport.diffusivity_used = D

    arch = architecture(model.land_sea_mask.land_fraction)
    land_fraction = on_architecture(CPU(), model.land_sea_mask.land_fraction)
    columns, weights, area = ring_diffusion_stencil(
        model.spectral_grid.grid, land_fraction, Float64(D), Float64(model.geometry.radius[])
    )

    heat_transport.max_rate = NF(maximum(sum(abs, weights, dims = 2)))
    heat_transport.columns = on_architecture(arch, columns)
    heat_transport.weights = on_architecture(arch, NF.(weights))
    heat_transport.area = on_architecture(arch, NF.(area))
    heat_transport.increment = on_architecture(arch, zeros(NF, length(area)))
    heat_transport.increment_next = on_architecture(arch, zeros(NF, length(area)))

    # explicit Euler is stable if Δt * (largest eigenvalue magnitude) < 2, which max_rate (Gershgorin) bounds,
    # substep such that the stability number is below max_stability, fixed here so that it is a constant later
    number = stability_number(heat_transport, model.time_stepping.Δt)
    heat_transport.substeps = max(1, ceil(Int, number / heat_transport.max_stability))
    return nothing
end

"""$(TYPEDSIGNATURES) Stability number of the explicit time stepping of `heat_transport` with time step `Δt` [s],
stable if below 1 (conservatively estimated, the actual limit is typically 1.5-2 times larger), before substepping."""
stability_number(heat_transport::OceanHeatDiffusion, Δt) = Δt * heat_transport.max_rate / 2

"""$(TYPEDSIGNATURES)
Precompute the conservative diffusion stencil on a ring `grid` with `land_fraction` and diffusivity `D` [m²/s] on a planet
of `radius` [m]. Returns `columns`, `weights` (both `npoints × n`, tendency `dT[i] = sum_m weights[i, m] * T[columns[i, m]]`)
and the cell `area` [m²]. Rows of `weights` sum to zero (uniform field stays uniform) and
`sum(area .* (1 .- land_fraction) .* dT) = 0` (conservative), for any `T`. For points without
any ocean the weights are zero."""
function ring_diffusion_stencil(grid::AbstractGrid, land_fraction, D::Real, radius::Real)
    Grid = RingGrids.nonparametric_type(typeof(grid))
    nlat_half = grid.nlat_half
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    rings = RingGrids.eachring(Grid, nlat_half)
    npoints = last(rings[end])
    londs, latds = RingGrids.get_londlatds(Grid, nlat_half)
    lats = [deg2rad(latds[first(ring)]) for ring in rings]          # latitude per ring
    @assert all(abs.(lats) .< π / 2 - 1.0e-9) "OceanHeatDiffusion does not support grids with rings at the poles"
    nlons = [length(ring) for ring in rings]
    lon0 = [deg2rad(londs[first(ring)]) for ring in rings]          # longitude of 1st point per ring
    lat_edges = [π / 2; (lats[1:(end - 1)] .+ lats[2:end]) ./ 2; -π / 2]

    ocean = clamp.(1 .- Float64.(land_fraction.data), 0, 1)   # ocean fraction
    rates = [Dict{Int, Float64}() for _ in 1:npoints]                # rates[a][b]: coefficient of T_b in dH_a/dt, in unit-sphere
    area = zeros(npoints)

    # heat exchange between cells a, b: rate_a += c (T_b - T_a), rate_b += c (T_a - T_b)
    add!(row, column, value) = rates[row][column] = get(rates[row], column, 0.0) + value
    function exchange!(a, b, c)
        (c == 0 || a == b) && return nothing
        for (i, j) in ((a, b), (b, a))
            rates[i][i] = get(rates[i], i, 0.0) - c
            rates[i][j] = get(rates[i], j, 0.0) + c
        end
        return nothing
    end

    for j in 1:nlat
        Δλ = 2π / nlons[j]
        area[rings[j]] .= Δλ * (sin(lat_edges[j]) - sin(lat_edges[j + 1])) * radius^2

        # zonal faces between neighbouring cells in the ring, face length Δφ, distance cos(φ) Δλ
        c0 = (lat_edges[j] - lat_edges[j + 1]) / (cos(lats[j]) * Δλ)
        for k in 1:nlons[j]
            a, b = first(rings[j]) + k - 1, first(rings[j]) + mod(k, nlons[j])
            exchange!(a, b, c0 * min(ocean[a], ocean[b]))
        end

        # meridional faces to the next ring south, interpolate both rings to the midpoints of the face segments
        j == nlat && continue
        faces = (rings[j], rings[j + 1])
        Δλs = (2π / nlons[j], 2π / nlons[j + 1])
        offsets = (lon0[j], lon0[j + 1])
        edges = Float64[]                                           # longitudes of cell edges, both rings
        for r in 1:2, k in 0:(nlons[j + r - 1] - 1)
            push!(edges, mod(offsets[r] + (k - 0.5) * Δλs[r], 2π))
        end
        sort!(edges)
        push!(edges, edges[1] + 2π)
        distance = lats[j] - lats[j + 1]
        for s in 1:(length(edges) - 1)
            length_s = edges[s + 1] - edges[s]
            length_s < 1.0e-12 && continue
            λ = (edges[s] + edges[s + 1]) / 2
            stencils = ntuple(2) do r                               # interpolation of ring r to λ, ocean-fraction weighted
                u = (λ - offsets[r]) / Δλs[r]
                k0 = floor(Int, u)
                frac = u - k0
                cells = (first(faces[r]) + mod(k0, nlons[j + r - 1]), first(faces[r]) + mod(k0 + 1, nlons[j + r - 1]))
                wgts = ((1 - frac) * ocean[cells[1]], frac * ocean[cells[2]])
                total = wgts[1] + wgts[2]
                inside = first(faces[r]) + mod(floor(Int, u + 0.5), nlons[j + r - 1])   # cell containing λ
                (cells, wgts ./ ifelse(total > 0, total, 1.0), inside)
            end
            w = min(ocean[stencils[1][3]], ocean[stencils[2][3]])
            c = w * cos(lat_edges[j + 1]) * length_s / distance
            c == 0 && continue
            # flux F = c (T_south - T_north) of interpolated temperatures, into the cell of the north ring and out of
            # the cell of the south ring that contain the segment
            north, south = stencils[1][3], stencils[2][3]
            for (cell, wgt) in zip(stencils[1][1], stencils[1][2])
                add!(north, cell, -c * wgt)
                add!(south, cell, c * wgt)
            end
            for (cell, wgt) in zip(stencils[2][1], stencils[2][2])
                add!(north, cell, c * wgt)
                add!(south, cell, -c * wgt)
            end
        end
    end

    # tendency of T = rate / (ocean fraction * area), with area on the unit sphere replaced by R² area, D / R²
    width = maximum(length.(rates))
    columns = repeat(collect(1:npoints), 1, width)
    weights = zeros(npoints, width)
    for i in 1:npoints
        ocean[i] > 0 || continue
        for (m, (col, val)) in enumerate(sort!(collect(rates[i])))
            columns[i, m] = col
            weights[i, m] = D * val / (ocean[i] * area[i])
        end
    end
    return columns, weights, area
end

function heat_transport!(dsst, vars, heat_transport::OceanHeatDiffusion, ocean, model)
    sst = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepper(model.time_stepping, :ocean), ocean)
    (; land_fraction) = model.land_sea_mask
    (; columns, weights, substeps) = heat_transport
    Δt = time_step(model.time_stepping, :ocean, vars.prognostic.clock)
    Δt_substep = Δt / substeps
    arch = architecture(dsst)

    # explicit substeps of the temperature increment, kept as increment to avoid cancellation of
    # the (small) increment against the temperature itself in low precision
    increment, increment_next = heat_transport.increment, heat_transport.increment_next
    increment .= 0
    for _ in 1:substeps
        launch!(
            arch, LinearWorkOrder, size(dsst), ocean_heat_diffusion_kernel!,
            increment_next, increment, sst, land_fraction, columns, weights, Δt_substep
        )
        increment, increment_next = increment_next, increment
    end
    launch!(arch, LinearWorkOrder, size(dsst), ocean_heat_tendency_kernel!, dsst, increment, land_fraction, inv(Δt))
    return nothing
end

@kernel inbounds = true function ocean_heat_diffusion_kernel!(increment_next, increment, sst, land_fraction, columns, weights, Δt)
    ij = @index(Global, Linear)         # every grid point ij
    if land_fraction[ij] < 1            # at least partially ocean
        tendency = zero(eltype(increment))
        for m in axes(weights, 2)
            column = columns[ij, m]
            tendency += weights[ij, m] * (sst[column] + increment[column])
        end
        increment_next[ij] = increment[ij] + Δt * tendency
    end
end

@kernel inbounds = true function ocean_heat_tendency_kernel!(dsst, increment, land_fraction, Δt⁻¹)
    ij = @index(Global, Linear)
    if land_fraction[ij] < 1
        dsst[ij] += Δt⁻¹ * increment[ij]
    end
end
