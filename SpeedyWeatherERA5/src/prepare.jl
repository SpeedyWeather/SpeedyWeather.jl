# Stage 2: ERA5 (native 0.25°, pressure levels) -> fields on the model's sigma layers and horizontal grid.
# Order of operations (surface pressure steps live in surface_pressure.jl):
#   1. surface pressure for the model (orography-consistent, or raw), see surface_pressure.jl
#   2. at native resolution, interpolate each column onto the model's sigma layers at p = σ·ps_target
#   3. truncate everything spectrally to the model resolution

const Γ_STD = 0.0065f0      # standard lapse rate [K/m], only used to extrapolate T below the lowest above-ground ERA5 level

"""
Interpolate `F` (nlon, nlat, nlev, on pressure levels `lev` [Pa], ascending) onto sigma layers `σ_half` of every column,
at p = σ·ps_target (the surface pressure the model layers sit on), linear in ln p. `ps_valid` is ERA5's own surface pressure,
levels with lev > ps_valid are below ERA5's ground. Each layer is sampled at its centre in σ. Levels below ground (lev > ps) are never used; if a target pressure lies below the lowest above-ground
level, T is extrapolated with the standard lapse rate and everything else held constant.
Returns the array on sigma layers and the fraction of samples that needed extrapolation."""
function pressure_to_sigma(F, lev, ps_valid, ps_target, σ_half; temperature = false, R = 287.0f0, g = 9.81f0)
    nlon, nlat, _ = size(F)
    nσ = length(σ_half) - 1
    out = zeros(Float32, nlon, nlat, nσ)
    nextrap = zeros(Int, nlat)
    κΓ = R * Γ_STD / g
    Threads.@threads for j in 1:nlat
        for i in 1:nlon
            p_s = ps_target[i, j]
            nv = max(searchsortedlast(lev, ps_valid[i, j]), 1)  # number of levels at or above ERA5's ground
            for k in 1:nσ
                σ = σ_half[k] + 0.5f0 * (σ_half[k + 1] - σ_half[k])     # layer centre
                pt = σ * p_s
                if pt <= lev[1]
                    val = F[i, j, 1]
                elseif pt <= lev[nv]
                    l = searchsortedlast(lev, pt)
                    if l >= nv
                        val = F[i, j, nv]
                    else
                        w = log(pt / lev[l]) / log(lev[l + 1] / lev[l])
                        val = (1 - w) * F[i, j, l] + w * F[i, j, l + 1]
                    end
                else
                    val = temperature ? F[i, j, nv] * (pt / lev[nv])^κΓ : F[i, j, nv]
                    nextrap[j] += 1
                end
                out[i, j, k] = val
            end
        end
    end
    return out, sum(nextrap) / (nlon * nlat * nσ)
end

"""Truncate a native-grid (nlon, nlat[, nlev]) array on the pole-less 0.25° grid spectrally to the model's resolution and
return it as a field on the model grid. `S_nat` is `SpectralTransform` of a `FullClenshawField` of that grid."""
function native_to_model(model, S_nat, a)
    sg = model.spectral_grid
    L = transform(FullClenshawField(Float32.(a), input_as = Matrix), S_nat)      # the native transform is planned in Float32
    Lm = ndims(a) == 2 ? zeros(sg.SpectralVariable2D, sg.spectrum) : zeros(sg.SpectralVariable3D, sg.spectrum, size(a, 3))
    copyto!(Lm, L)                                              # copies the overlapping (= truncated) set of coefficients
    return transform(Lm, model.spectral_transform)
end

"""
Ocean-only box average of a native (nlon, nlat) field on the pole-less 0.25° grid onto the model's FullGaussian grid:
every valid (non-NaN) native cell contributes to the model cell it falls in, weighted by cos(lat); no spectral ringing at coasts
or ice edges. Model cells without any valid native cell are filled with the mean of their filled neighbours (`fill_value`
= `nothing`) or with `fill_value`."""
function masked_to_model(model, a; fill_value = nothing)
    nlh = model.spectral_grid.nlat_half
    nlon, nlat = 4nlh, 2nlh
    latm = RingGrids.get_latd(FullGaussianGrid, nlh)
    latn = RingGrids.get_latd(FullClenshawGrid, 360)
    edges = [90.0; (latm[1:(end - 1)] .+ latm[2:end]) ./ 2; -90.0]                  # descending model cell edges
    jrow = [clamp(findfirst(<(l), edges) - 1, 1, nlat) for l in latn]
    icol = [mod(round(Int, (i - 1) * 0.25 / (360 / nlon)), nlon) + 1 for i in 1:size(a, 1)]
    num, den = zeros(Float64, nlon, nlat), zeros(Float64, nlon, nlat)
    for jn in axes(a, 2), i in axes(a, 1)
        x = a[i, jn]
        isnan(x) && continue
        w = cosd(latn[jn])
        num[icol[i], jrow[jn]] += w * x
        den[icol[i], jrow[jn]] += w
    end
    out = num ./ den                                            # NaN where den == 0
    if fill_value === nothing
        while any(isnan, out)                                   # grow into empty cells from valid neighbours
            new = copy(out)
            for j in 1:nlat, i in 1:nlon
                isnan(out[i, j]) || continue
                nb = [out[mod1(i + di, nlon), j + dj] for di in -1:1, dj in -1:1 if (di != 0 || dj != 0) && 1 <= j + dj <= nlat]
                v = filter(!isnan, nb)
                isempty(v) || (new[i, j] = sum(v) / length(v))
            end
            out = new
        end
    else
        out[isnan.(out)] .= fill_value
    end
    return FullGaussianField(Float32.(out), input_as = Matrix)
end

"""
Prepare the atmospheric fields of the cached ERA5 file `path` for `model` (a PrimitiveWetModel): returns
(; temperature, humidity, u, v, pressure, sst, sea_ice, info) with fields on the model grid, `pressure` = log(ps [Pa]) consistent with
the MODEL's orography. Pole rows of ERA5 are dropped so that the 0.25° data is exactly a FullClenshawGrid(360)."""
function prepare_era5(model, path; target = :model, surface_pressure = :consistent)
    sg = model.spectral_grid
    S = model.spectral_transform
    g = model.planet.gravity
    R = model.atmosphere.R_dry
    σ_half = Array(model.geometry.σ_levels_half)

    ds = NCDataset(path)
    lev = Float32.(ds["level"][:]) .* 100                       # hPa -> Pa, ascending
    rd2(v) = Float32.(ds[v][:, 2:720])                          # drop pole rows
    rd3(v) = Float32.(ds[v][:, 2:720, :])
    ps, zs = rd2("surface_pressure"), rd2("geopotential_at_surface")
    sst_native, ice_native = rd2("sea_surface_temperature"), rd2("sea_ice_cover")
    T, q, u, v = rd3("temperature"), rd3("specific_humidity"), rd3("u_component_of_wind"), rd3("v_component_of_wind")
    close(ds)

    # native spectral transform, planned once; truncate to the model's spectrum, transform onto the model grid
    S_nat = SpectralTransform(FullClenshawField(zs, input_as = Matrix))
    on_model(a) = native_to_model(model, S_nat, a)

    # surface pressure: orography-consistent (default, src/surface_pressure.jl) or ERA5's raw ps truncated (naive, for comparison only)
    pressure = surface_pressure == :raw ? raw_pressure(on_model, ps) :
               orography_consistent_pressure(model, on_model, T, lev, ps, zs; R, g)

    # the pressure the model's layers sit on, as a smooth field on the native grid (model-resolution ps expanded onto 0.25°)
    if target == :model
        Lps = transform(pressure, S)
        Lnat = zeros(ComplexF32, S_nat.spectrum)               # native spectrum, zero beyond the model truncation
        copyto!(Lnat, Lps)
        ps_target = exp.(Float32.(Matrix(transform(Lnat, S_nat))))
    else                                                        # :era5, sample on ERA5's own terrain-following surfaces
        ps_target = ps
    end

    # vertical: native columns -> sigma layers at p = σ·ps_target
    Tσ, fT = pressure_to_sigma(T, lev, ps, ps_target, σ_half; temperature = true, R, g)
    qσ, _ = pressure_to_sigma(q, lev, ps, ps_target, σ_half)
    uσ, _ = pressure_to_sigma(u, lev, ps, ps_target, σ_half)
    vσ, fuv = pressure_to_sigma(v, lev, ps, ps_target, σ_half)

    temperature, humidity = on_model(Tσ), on_model(qσ)
    humidity .= max.(humidity, 0)
    sst = masked_to_model(model, sst_native)
    sea_ice = masked_to_model(model, ice_native; fill_value = 0.0f0)
    info = (; extrapolated_fraction_T = fT, extrapolated_fraction_uv = fuv, target, surface_pressure)
    return (; temperature, humidity, u = on_model(uσ), v = on_model(vσ), pressure, sst, sea_ice, info)
end
