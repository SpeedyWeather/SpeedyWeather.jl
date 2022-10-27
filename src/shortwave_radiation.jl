"""
    function shortwave_radiation!(
        column::ColumnVariables{NF}, model::PrimitiveEquationModel
    )

Compute air temperature tendencies from shortwave radiation for an atmospheric column.
For more details see http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf
"""
function shortwave_radiation!(
    column::ColumnVariables{NF}, model::PrimitiveEquationModel
) where {NF<:AbstractFloat}
    @unpack humid, sat_vap_pres, dry_static_energy, geopot, norm_pres = column
    @unpack cp, p0, = model.parameters
    @unpack σ_levels_thick = model.geometry
    @unpack gravity = model.constants

    sol_oz!(column, model)

    column.rel_hum = humid ./ sat_vap_pres
    column.grad_dry_static_energy =
        (dry_static_energy[end - 1] - dry_static_energy[end]) ./
        (geopot[end - 1] - geopot[end])
    cloud!(column, model)

    radsw!(column, model)

    prog_sp_inv = 1.0 / norm_pres
    grdsig = gravity ./ (σ_levels_thick .* p0)

    for k in eachlayer(column)
        column.temp_tend[k] += column.tend_t_rsw[k] * prog_sp_inv * grdsig[k] / cp
    end
end

"""
    function solar!(column::ColumnVariables{NF})

Compute average daily flux of solar radiation for an atmospheric column,
from Hartmann (1994).
"""
function solar!(column::ColumnVariables{NF}) where {NF<:AbstractFloat}
    @unpack tyear, csol, lat = column

    # Compute cosine and sine of latitude
    clat = cos(lat * π / 180.0)
    slat = sin(lat * π / 180.0)

    # 1.0 Compute declination angle and Earth - Sun distance factor
    pigr = 2.0 * asin(1.0)
    α = 2.0 * pigr * tyear

    ca1 = cos(α)
    sa1 = sin(α)
    ca2 = ca1 * ca1 - sa1 * sa1
    sa2 = 2.0 * sa1 * ca1
    ca3 = ca1 * ca2 - sa1 * sa2
    sa3 = sa1 * ca2 + sa2 * ca1

    decl =
        0.006918 - 0.399912 * ca1 + 0.070257 * sa1 - 0.006758 * ca2 + 0.000907 * sa2 -
        0.002697 * ca3 + 0.001480 * sa3

    fdis = 1.000110 + 0.034221 * ca1 + 0.001280 * sa1 + 0.000719 * ca2 + 0.000077 * sa2

    cdecl = cos(decl)
    sdecl = sin(decl)
    tdecl = sdecl / cdecl

    # 2.0 Compute daily - average insolation at the atm. top
    csolp = csol / pigr

    ch0 = min(1.0, max(-1.0, -tdecl * slat / clat))
    h0 = acos(ch0)
    sh0 = sin(h0)

    column.topsr = csolp * fdis * (h0 * slat * sdecl + sh0 * clat * cdecl)
end

"""
    function sol_oz!(
        column::ColumnVariables{NF}, model::PrimitiveEquationModel
    )

Compute solar radiation parametres for an atmospheric column.
"""
function sol_oz!(
    column::ColumnVariables{NF}, model::PrimitiveEquationModel
) where {NF<:AbstractFloat}
    @unpack tyear, lat = column
    @unpack solc, epssw = model.parameters

    # Compute cosine and sine of latitude
    clat = cos(lat * π / 180.0)
    slat = sin(lat * π / 180.0)

    # α = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
    α = 4.0 * asin(1.0) * (tyear + 10.0 / 365.0)
    dα = 0.0

    coz1 = 1.0 * max(0.0, cos(α - dα))
    coz2 = 1.8
    azen = 1.0
    nzen = 2

    rzen = -cos(α) * 23.45 * asin(1.0) / 90.0
    czen = cos(rzen)
    szen = sin(rzen)

    fs0 = 6.0

    # Solar radiation at the top
    column.csol = 4.0 * solc
    column.topsr = solar!(column)
    column.fsol = column.topsr

    flat2 = 1.5 * slat^2 - 0.5

    # Ozone depth in upper and lower stratosphere
    column.ozupp = 0.5 * epssw
    column.ozone = 0.4 * epssw * (1.0 + coz1 * slat + coz2 * flat2)

    # Zenith angle correction to (downward) absorptivity
    column.zenit = 1.0 + azen * (1.0 - (clat * czen + slat * szen))^nzen

    # Ozone absorption in upper and lower stratosphere
    column.ozupp = column.fsol * column.ozupp * column.zenit
    column.ozone = column.fsol * column.ozone * column.zenit

    # Polar night cooling in the stratosphere
    column.stratz = max(fs0 - column.fsol, 0.0)
end

"""
    function cloud!(
        column::ColumnVariables{NF}, model::PrimitiveEquationModel
    )

Compute shortwave radiation cloud contibutions for an atmospheric column.
"""
function cloud!(
    column::ColumnVariables{NF}, model::PrimitiveEquationModel
) where {NF<:AbstractFloat}
    @unpack rhcl1, rhcl2, rrcl, qcl, pmaxcl, wpcl, gse_s1, gse_s0, clsmax, clsminl =
        model.parameters
    @unpack humid, rel_hum, grad_dry_static_energy, precip_convection, precip_large_scale,
        cloud_top, nlev, fmask = column

    # 1.0 Cloud cover, defined as the sum of:
    # - a term proportional to the square - root of precip. rate
    # - a quadratic function of the max. relative humidity
    #    in tropospheric layers above pbl where q > prog_qcl :
    #     ( = 0 for rhmax < rhcl1, = 1 for rhmax > rhcl2 )
    #    Cloud - top level: defined as the highest (i.e. least sigma)
    #      between the top of convection / condensation and
    #      the level of maximum relative humidity.
    if (rel_hum[nlev - 1] > rhcl1)
        column.cloudc = rel_hum[nlev - 1] - rhcl1
        column.icltop = nlev - 1
    else
        column.cloudc = 0.0
        column.icltop = nlev + 1
    end

    for k in 3:(nlev - 2)
        drh = rel_hum[k] - rhcl1
        if (drh > column.cloudc) & (humid[k] > qcl)
            column.cloudc = drh
            column.icltop = k
        end
    end

    pr1 = min(pmaxcl, 86.4 * (precip_convection + precip_large_scale))
    column.cloudc = min(1.0, wpcl * sqrt(pr1) + min(1.0, column.cloudc * rrcl)^2.0)
    column.icltop = min(cloud_top, column.icltop)

    # 2.0  Equivalent specific humidity of clouds
    column.qcloud = humid[nlev - 1]

    # 3.0 Stratiform clouds at the top of pbl
    clfact = 1.2
    rgse = 1.0 / (gse_s1 - gse_s0)

    # Stratocumulus clouds over sea
    fstab = max(0.0, min(1.0, rgse * (grad_dry_static_energy - gse_s0)))
    column.clstr = fstab * max(clsmax - clfact * column.cloudc, 0.0)

    # Stratocumulus clouds over land
    clstrl = max(column.clstr, clsminl) * rel_hum[nlev]
    column.clstr = column.clstr + fmask * (clstrl - column.clstr)
end

"""
    function radsw!(
        column::ColumnVariables{NF}, model::PrimitiveEquationModel
    )

Compute shortwave radiation fluxes for an atmospheric column.
"""
function radsw!(
    column::ColumnVariables{NF}, model::PrimitiveEquationModel
) where {NF<:AbstractFloat}
    @unpack norm_pres,
    humid, icltop, cloudc, clstr, ozupp, ozone, zenit, stratz, fsol, qcloud, albsfc,
        nlev = column
    @unpack σ_levels_full, σ_levels_thick = model.geometry

    @unpack albcl, albcls, abscl1, abscl2, absdry, absaer, abswv1, abswv2, ablwin, ablco2, ablwv1,
        ablwv2, ablcl2, ablcl1, epslw = model.parameters

    # Locals variables
    flux::Vector{NF} = fill(NaN, 2)
    stratc::Vector{NF} = fill(NaN, 2)

    nl1 = nlev - 1

    fband2 = 0.05
    fband1 = 1.0 - fband2

    # 1.0 Initialization
    column.tau2 = fill!(column.tau2, 0)

    # Change to ensure only ICLTOP < = NLEV used
    if (icltop <= nlev)
        column.tau2[icltop, 3] = albcl * cloudc
    end
    column.tau2[nlev, 3] = albcls * clstr

    # 2.0 Shortwave transmissivity:
    # function of layer mass, ozone (in the statosphere),
    # abs. humidity and cloud cover (in the troposphere)
    psaz = norm_pres * zenit
    acloud = cloudc * min(abscl1 * qcloud, abscl2)

    deltap = psaz * σ_levels_thick[1]
    column.tau2[1, 1] = exp(-deltap * absdry)

    for k in 2:nl1
        abs1 = absdry + absaer * σ_levels_full[k]^2
        deltap = psaz * σ_levels_thick[k]
        if (k >= icltop)
            column.tau2[k, 1] = exp(-deltap * (abs1 + abswv1 * humid[k] + acloud))
        else
            column.tau2[k, 1] = exp(-deltap * (abs1 + abswv1 * humid[k]))
        end
    end

    abs1 = absdry + absaer * σ_levels_full[nlev]^2
    deltap = psaz * σ_levels_thick[nlev]
    column.tau2[nlev, 1] = exp(-deltap * (abs1 + abswv1 * humid[nlev]))

    for k in 2:nlev
        deltap = psaz * σ_levels_thick[k]
        column.tau2[k, 2] = exp(-deltap * abswv2 * humid[k])
    end

    # 3.0 Shortwave downward flux
    # 3.1 Initialization of fluxes
    column.tsr = fsol
    flux[1] = fsol * fband1
    flux[2] = fsol * fband2

    # 3.2 Ozone and dry - air absorption in the stratosphere
    k = 1
    column.tend_t_rsw[k] = flux[1]
    flux[1] = column.tau2[k, 1] * (flux[1] - ozupp * norm_pres)
    column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]

    k = 2
    column.tend_t_rsw[k] = flux[1]
    flux[1] = column.tau2[k, 1] * (flux[1] - ozone * norm_pres)
    column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]

    # 3.3  Absorption and reflection in the troposphere
    for k in 3:nlev
        column.tau2[k, 3] = flux[1] * column.tau2[k, 3]
        flux[1] = flux[1] - column.tau2[k, 3]
        column.tend_t_rsw[k] = flux[1]
        flux[1] = column.tau2[k, 1] * flux[1]
        column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]
    end

    for k in 2:nlev
        column.tend_t_rsw[k] = column.tend_t_rsw[k] + flux[2]
        flux[2] = column.tau2[k, 2] * flux[2]
        column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[2]
    end

    # 4.0 Shortwave upward flux
    # 4.1  Absorption and reflection at the surface
    column.ssrd = flux[1] + flux[2]
    flux[1] = flux[1] * albsfc
    column.ssr = column.ssrd - flux[1]

    # 4.2  Absorption of upward flux
    for k in nlev:-1:1
        column.tend_t_rsw[k] = column.tend_t_rsw[k] + flux[1]
        flux[1] = column.tau2[k, 1] * flux[1]
        column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]
        flux[1] = flux[1] + column.tau2[k, 3]
    end

    # 4.3  Net solar radiation = incoming - outgoing
    column.tsr = column.tsr - flux[1]

    # 5.0  Initialization of longwave radiation model
    # 5.1  Longwave transmissivity:
    #      function of layer mass, abs. humidity and cloud cover.

    #    Cloud - free levels (stratosphere + PBL)
    k = 1
    deltap = norm_pres * σ_levels_thick[k]
    column.tau2[k, 1] = exp(-deltap * ablwin)
    column.tau2[k, 2] = exp(-deltap * ablco2)
    column.tau2[k, 3] = 1.0
    column.tau2[k, 4] = 1.0

    for k in 2:(nlev - 2):nlev
        deltap = norm_pres * σ_levels_thick[k]
        column.tau2[k, 1] = exp(-deltap * ablwin)
        column.tau2[k, 2] = exp(-deltap * ablco2)
        column.tau2[k, 3] = exp(-deltap * ablwv1 * humid[k])
        column.tau2[k, 4] = exp(-deltap * ablwv2 * humid[k])
    end

    # Cloudy layers (free troposphere)
    acloud = cloudc * ablcl2

    for k in 3:nl1
        deltap = norm_pres * σ_levels_thick[k]
        if (k < icltop)
            acloud1 = acloud
        else
            acloud1 = ablcl1 * cloudc
        end
        column.tau2[k, 1] = exp(-deltap * (ablwin + acloud1))
        column.tau2[k, 2] = exp(-deltap * ablco2)
        column.tau2[k, 3] = exp(-deltap * max(ablwv1 * humid[k], acloud))
        column.tau2[k, 4] = exp(-deltap * max(ablwv2 * humid[k], acloud))
    end

    # 5.2  Stratospheric correction terms
    eps1 = epslw / (σ_levels_thick[1] + σ_levels_thick[2])
    stratc[1] = stratz * norm_pres
    stratc[2] = eps1 * norm_pres
end
