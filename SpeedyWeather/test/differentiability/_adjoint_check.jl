# Develop + FD-validate the CHUNKED hand-written transform! adjoints (NO Enzyme → compiles fast).
# These functions mirror what goes into SpeedyTransformsEnzymeExt. Validated on a batched S
# (1 chunk) and a chunked S (transform_batch=[1] → per-layer chunks) against FiniteDifferences.
# Run: julia --project=SpeedyWeather/test/differentiability <file>
using SpeedyWeather
using FiniteDifferences
import SpeedyWeather.SpeedyTransforms:
    SpectralTransform, transform!, _fourier!, _legendre!, wrapped_view,
    _largest_planned_batch, ColumnScratchMemory
import Random

println("Julia ", VERSION)
Random.seed!(7)

# adjoint_scale / rfft_adjoint_scale copied from SpeedyTransformsEnzymeExt (only defined there)
function rfft_adjoint_scale(n_freq::Int, n_real::Int)
    if iseven(n_real)
        return [1 < i < n_freq ? 2 : 1 for i in 1:n_freq]
    else
        return [1 < i ? 2 : 1 for i in 1:n_freq]
    end
end
function adjoint_scale(S::SpectralTransform)
    (; nlons) = S
    (; nlat_half) = S.grid
    rfft_plans_1D = S.rfft_plans[1]
    nfreqs = [rfft_plan.osz[1] for rfft_plan in rfft_plans_1D]
    scale = zeros(Int, maximum(nfreqs), 1, nlat_half)
    for i in 1:nlat_half
        scale[1:nfreqs[i], 1, i] = rfft_adjoint_scale(nfreqs[i], nlons[i])
    end
    return scale
end

# ---- CHUNKED adjoint of spec->grid transform! w.r.t coeffs (accumulate into coeffs_bar) ----
function spec2grid_pullback!(coeffs_bar, field_bar, S; unscale_coslat::Bool = false)
    NF = eltype(S)
    (; nlat_half) = S.grid
    K = size(field_bar, 2)
    Kb = _largest_planned_batch(K, S)
    scale = adjoint_scale(S)
    dOmega = reshape(view(S.solid_angles, 1:nlat_half), 1, 1, :)
    clat = reshape(view(S.coslat⁻¹, 1:nlat_half), 1, 1, :)
    c = 1
    while c <= K
        ce = min(c + Kb - 1, K)
        chunk = c:ce
        Kc = ce - c + 1
        dg_n = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        dg_s = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        _fourier!(dg_n, dg_s, wrapped_view(field_bar, :, chunk), S)   # forward FFT
        dg_n .*= scale
        dg_s .*= scale
        if unscale_coslat
            dg_n .*= clat
            dg_s .*= clat
        end
        dg_n ./= dOmega
        dg_s ./= dOmega
        col = ColumnScratchMemory(zeros(Complex{NF}, Kc), zeros(Complex{NF}, Kc))
        cbar = zeros(Complex{NF}, S.spectrum, Kc)
        _legendre!(cbar, dg_n, dg_s, col, S)                          # forward Legendre
        wrapped_view(coeffs_bar, :, chunk).data .+= cbar.data
        c = ce + 1
    end
    return coeffs_bar
end

# ---- CHUNKED adjoint of grid->spec transform! w.r.t field (accumulate into field_bar) ----
function grid2spec_pullback!(field_bar, coeffs_bar, S)
    NF = eltype(S)
    (; nlat_half) = S.grid
    K = size(coeffs_bar, 2)
    Kb = _largest_planned_batch(K, S)
    scale = adjoint_scale(S)
    dOmega = reshape(view(S.solid_angles, 1:nlat_half), 1, 1, :)
    c = 1
    while c <= K
        ce = min(c + Kb - 1, K)
        chunk = c:ce
        Kc = ce - c + 1
        df_n = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        df_s = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        col = ColumnScratchMemory(zeros(Complex{NF}, Kc), zeros(Complex{NF}, Kc))
        _legendre!(df_n, df_s, wrapped_view(coeffs_bar, :, chunk), col, S)   # inverse Legendre
        df_n .*= dOmega
        df_s .*= dOmega
        df_n ./= scale
        df_s ./= scale
        fbar = zeros(NF, S.grid, Kc)
        _fourier!(fbar, df_n, df_s, S)                                # inverse FFT
        wrapped_view(field_bar, :, chunk).data .+= fbar.data
        c = ce + 1
    end
    return field_bar
end

# ================= FD validation =================
NL = 4
spectral_grid = SpectralGrid(trunc = 5, nlayers = NL)
spectrum = spectral_grid.spectrum
grid = spectral_grid.grid

coeffs0 = rand(ComplexF32, spectrum, NL)
field0 = rand(Float32, grid, NL)
field_bar0 = rand(Float32, grid, NL)
coeffs_bar0 = rand(ComplexF32, spectrum, NL)

for (name, batch) in (("batched [1,4]", [1, NL]), ("chunked [1]", [1]))
    S = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = batch)

    # spec->grid, vjp w.r.t coeffs
    primal_s2g(cvec) = begin
        coeffs = deepcopy(coeffs0)
        coeffs.data .= reshape(reinterpret(ComplexF32, copy(cvec)), size(coeffs.data))
        field = zeros(Float32, grid, NL)
        transform!(field, coeffs, deepcopy(S.scratch_memory), S)
        vec(copy(field.data))
    end
    cvec0 = collect(reinterpret(Float32, vec(copy(coeffs0.data))))
    fd = FiniteDifferences.j′vp(central_fdm(5, 1), primal_s2g, vec(field_bar0.data), cvec0)[1]
    fd_c = reshape(reinterpret(ComplexF32, fd), size(coeffs0.data))
    mine = spec2grid_pullback!(zeros(ComplexF32, spectrum, NL), field_bar0, S)
    println(rpad(name, 14), " spec->grid rel = ",
        maximum(abs, mine.data .- fd_c) / maximum(abs, fd_c))

    # grid->spec, vjp w.r.t field
    primal_g2s(fvec) = begin
        field = zeros(Float32, grid, NL)
        field.data .= reshape(copy(fvec), size(field.data))
        coeffs = zeros(ComplexF32, spectrum, NL)
        transform!(coeffs, field, deepcopy(S.scratch_memory), S)
        collect(reinterpret(Float32, vec(copy(coeffs.data))))
    end
    fvec0 = vec(copy(field0.data))
    cbarvec = collect(reinterpret(Float32, vec(copy(coeffs_bar0.data))))
    fd2 = FiniteDifferences.j′vp(central_fdm(5, 1), primal_g2s, cbarvec, fvec0)[1]
    fd_f = reshape(fd2, size(field0.data))
    mine2 = grid2spec_pullback!(zeros(Float32, grid, NL), coeffs_bar0, S)
    println(rpad(name, 14), " grid->spec rel = ",
        maximum(abs, mine2.data .- fd_f) / maximum(abs, fd_f))
end
