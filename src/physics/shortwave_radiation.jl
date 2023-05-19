"""
Parameters for radiation parameterizations.
"""
@kwdef struct RadiationCoefs{NF<:Real} <: Coefficients
    epslw::NF = 0.05    # Fraction of blackbody spectrum absorbed/emitted by PBL only
    emisfc::NF = 0.98   # Longwave surface emissivity

    p0::Real = 1e5      # Reference pressure (Pa)

    # Shortwave radiation: sol_oz
    solc::NF = 342.0    # Solar constant (area averaged) [W/m^2]
    epssw::NF = 0.020   # Fraction of incoming solar radiation absorbed by ozone

    # Shortwave radiation: cloud
    rhcl1::NF = 0.30    # Relative humidity threshold corresponding to cloud cover = 0
    rhcl2::NF = 1.00    # Relative humidity correponding to cloud cover = 1
    rrcl::NF = 1 / (rhcl2 - rhcl1)
    qcl::NF = 0.20      # Specific humidity threshold for cloud cover
    pmaxcl::NF = 10.0   # Maximum value of precipitation (mm/day) contributing to cloud cover
    wpcl::NF = 0.2      # Cloud cover weight for the square-root of precipitation (for p = 1 mm/day)
    gse_s1::NF = 0.40   # Gradient of dry static energy corresponding to stratiform cloud cover = 1
    gse_s0::NF = 0.25   # Gradient of dry static energy corresponding to stratiform cloud cover = 0
    clsmax::NF = 0.60   # Maximum stratiform cloud cover
    clsminl::NF = 0.15  # Minimum stratiform cloud cover over land (for RH = 1)

    # Shortwave radiation: radsw
    albcl::NF = 0.43    # Cloud albedo (for cloud cover = 1)
    albcls::NF = 0.50   # Stratiform cloud albedo (for st. cloud cover = 1)
    abscl1::NF = 0.015  # Absorptivity of clouds (visible band, maximum value)
    abscl2::NF = 0.15   # Absorptivity of clouds (visible band, for dq_base = 1 g/kg)

    absdry::NF = 0.033  # Absorptivity of dry air (visible band)
    absaer::NF = 0.033  # Absorptivity of aerosols (visible band)
    abswv1::NF = 0.022  # Absorptivity of water vapour (visible band, for dq = 1 g/kg)
    abswv2::NF = 15.0   # Absorptivity of water vapour (near IR band, for dq = 1 g/kg)

    ablwin::NF = 0.3    # Absorptivity of air in "window" band
    ablco2::NF = 6.0    # Absorptivity of air in CO2 band
    ablwv1::NF = 0.7    # Absorptivity of water vapour in H2O band 1 (weak), (for dq = 1 g/kg)
    ablwv2::NF = 50.0   # Absorptivity of water vapour in H2O band 2 (strong), (for dq = 1 g/kg)
    ablcl2::NF = 0.6    # Absorptivity of "thin" upper clouds in window and H2O bands
    ablcl1::NF = 12.0   # Absorptivity of "thick" clouds in window band (below cloud top)
end 

"""
    function shortwave_radiation!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    )

Compute air temperature tendencies from shortwave radiation for an atmospheric column.
For more details see http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf
"""
function shortwave_radiation!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}
    (; humid, sat_vap_pres, dry_static_energy, geopot, norm_pres ) = column
    (; cₚ ) = model.parameters
    (; p0 ) = model.parameters.radiation_coefs
    (; σ_levels_thick ) = model.geometry
    (; gravity ) = model.constants

    sol_oz!(column, model)

    column.rel_hum .= humid ./ sat_vap_pres
    column.grad_dry_static_energy = 
        (dry_static_energy[end - 1] - dry_static_energy[end]) /
        (geopot[end - 1] - geopot[end])
    cloud!(column, model)

    radsw!(column, model)

    for k in eachlayer(column)
        column.temp_tend[k] += 
            column.tend_t_rsw[k] * inv(norm_pres) *
            (gravity / (σ_levels_thick[k] * p0)) /cₚ
    end
end

"""
    function solar!(column::ColumnVariables{NF})

Compute average daily flux of solar radiation for an atmospheric column,
from Hartmann (1994).
"""
function solar!(column::ColumnVariables{NF}) where {NF<:AbstractFloat}
    (; tyear, csol, latd ) = column

    # Compute cosine and sine of latitude
    clat = cos(latd)
    slat = sin(latd)

    # 1.0 Compute declination angle and Earth - Sun distance factor
    pigr = 2 * asin(1)
    α = 2 * pigr * tyear

    ca1 = cos(α)
    sa1 = sin(α)
    ca2 = ca1 * ca1 - sa1 * sa1
    sa2 = 2 * sa1 * ca1
    ca3 = ca1 * ca2 - sa1 * sa2
    sa3 = sa1 * ca2 + sa2 * ca1

    decl =
        NF(0.006918) - NF(0.399912) * ca1 + NF(0.070257) * sa1 - NF(0.006758) * ca2 + NF(0.000907) * sa2 -
        NF(0.002697) * ca3 + NF(0.001480) * sa3

    fdis = NF(1.000110) + NF(0.034221) * ca1 + NF(0.001280) * sa1 + NF(0.000719) * ca2 + NF(0.000077) * sa2

    cdecl = cos(decl)
    sdecl = sin(decl)
    tdecl = sdecl / cdecl

    # 2.0 Compute daily - average insolation at the atm. top
    csolp = csol / pigr

    ch0 = min(1, max(-1, -tdecl * slat / clat))
    h0 = acos(ch0)
    sh0 = sin(h0)

    column.topsr = csolp * fdis * (h0 * slat * sdecl + sh0 * clat * cdecl)
end

"""
    function sol_oz!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    )

Compute solar radiation parametres for an atmospheric column.
"""
function sol_oz!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}
    (; tyear, latd ) = column
    (; solc, epssw ) = model.parameters.radiation_coefs

    # Compute cosine and sine of latitude
    clat = cos(latd)
    slat = sin(latd)

    # α = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
    α = 4 * asin(1) * (tyear + 10 / 365)
    dα = NF(0.)

    coz1 = 1 * max(0, cos(α - dα))
    coz2 = NF(1.8)
    azen = 1
    nzen = 2

    rzen = -cos(α) * NF(23.45) * asin(1.0) / 90
    czen = cos(rzen)
    szen = sin(rzen)

    fs0 = 6

    # Solar radiation at the top
    column.csol = 4 * solc
    column.topsr = solar!(column)
    column.fsol = column.topsr

    flat2 = NF(1.5) * slat^2 - NF(0.5)

    # Ozone depth in upper and lower stratosphere
    column.ozupp = NF(0.5) * epssw
    column.ozone = NF(0.4) * epssw * (1 + coz1 * slat + coz2 * flat2)

    # Zenith angle correction to (downward) absorptivity
    column.zenit = 1 + azen * (1 - (clat * czen + slat * szen))^nzen

    # Ozone absorption in upper and lower stratosphere
    column.ozupp = column.fsol * column.ozupp * column.zenit
    column.ozone = column.fsol * column.ozone * column.zenit

    # Polar night cooling in the stratosphere
    column.stratz = max(fs0 - column.fsol, 0)
end

"""
    function cloud!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    )

Compute shortwave radiation cloud contibutions for an atmospheric column.
"""
function cloud!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}
    (; rhcl1, rhcl2, rrcl, qcl, pmaxcl ) = model.parameters.radiation_coefs
    (; wpcl, gse_s1, gse_s0, clsmax, clsminl ) = model.parameters.radiation_coefs
    (; humid, rel_hum, grad_dry_static_energy, precip_convection ) = column
    (; precip_large_scale, cloud_top, nlev, fmask ) = column
    (; n_stratosphere_levels ) = model.geometry

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
        column.cloudc = NF(0)
        column.icltop = nlev + 1
    end

    for k in n_stratosphere_levels+1:(nlev - n_stratosphere_levels)
        drh = rel_hum[k] - rhcl1
        if (drh > column.cloudc) & (humid[k] > qcl)
            column.cloudc = drh
            column.icltop = k
        end
    end

    pr1 = min(pmaxcl, NF(86.4) * (precip_convection + precip_large_scale))
    column.cloudc = min(1, wpcl * sqrt(pr1) + min(1, column.cloudc * rrcl)^2)
    column.icltop = min(cloud_top, column.icltop)

    # 2.0  Equivalent specific humidity of clouds
    column.qcloud = humid[nlev - 1]

    # 3.0 Stratiform clouds at the top of pbl
    clfact = NF(1.2)
    rgse = 1 / (gse_s1 - gse_s0)

    # Stratocumulus clouds over sea
    fstab = max(0, min(1, rgse * (grad_dry_static_energy - gse_s0)))
    column.clstr = fstab * max(clsmax - clfact * column.cloudc, 0)

    # Stratocumulus clouds over land
    clstrl = max(column.clstr, clsminl) * rel_hum[nlev]
    column.clstr = column.clstr + fmask * (clstrl - column.clstr)
end

"""
    function radsw!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    )

Compute shortwave radiation fluxes for an atmospheric column.
"""
function radsw!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}
    (; norm_pres, humid, icltop, cloudc, clstr, ozupp, ozone ) = column
    (; zenit, stratz, fsol, qcloud, albsfc, nlev ) = column
    (; σ_levels_full, σ_levels_thick, n_stratosphere_levels ) = model.geometry
    (; albcl, albcls, abscl1, abscl2, absdry, absaer ) = model.parameters.radiation_coefs
    (; abswv1, abswv2, ablwin, ablco2, ablwv1 ) = model.parameters.radiation_coefs
    (; ablwv2, ablcl2, ablcl1, epslw ) = model.parameters.radiation_coefs

    # Locals variables
    sbands_flux = 2
    flux::Vector{NF} = fill(NaN, sbands_flux)
    stratc::Vector{NF} = fill(NaN, n_stratosphere_levels)

    nl1 = nlev - 1

    fband2 = NF(0.05)
    fband1 = 1 - fband2

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

    for k in n_stratosphere_levels:nl1
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

    for k in n_stratosphere_levels:nlev
        deltap = psaz * σ_levels_thick[k]
        column.tau2[k, 2] = exp(-deltap * abswv2 * humid[k])
    end

    # 3.0 Shortwave downward flux
    # 3.1 Initialization of fluxes
    column.tsr = fsol
    flux[1] = fsol * fband1
    flux[2] = fsol * fband2

    # 3.2 Ozone and dry - air absorption in the stratosphere
    for k in 1:n_stratosphere_levels
        if k == 1 ozone_tmp = ozupp else ozone_tmp = ozone end
        column.tend_t_rsw[k] = flux[1]
        flux[1] = column.tau2[k, 1] * (flux[1] - ozone_tmp * norm_pres)
        column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]
    end
    
    # 3.3  Absorption and reflection in the troposphere
    for k in n_stratosphere_levels+1:nlev
        column.tau2[k, 3] = flux[1] * column.tau2[k, 3]
        flux[1] = flux[1] - column.tau2[k, 3]
        column.tend_t_rsw[k] = flux[1]
        flux[1] = column.tau2[k, 1] * flux[1]
        column.tend_t_rsw[k] = column.tend_t_rsw[k] - flux[1]
    end

    for k in n_stratosphere_levels:nlev
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
    column.tau2[k, 3] = 1.
    column.tau2[k, 4] = 1.

    for k in n_stratosphere_levels:(nlev - n_stratosphere_levels):nlev
        deltap = norm_pres * σ_levels_thick[k]
        column.tau2[k, 1] = exp(-deltap * ablwin)
        column.tau2[k, 2] = exp(-deltap * ablco2)
        column.tau2[k, 3] = exp(-deltap * ablwv1 * humid[k])
        column.tau2[k, 4] = exp(-deltap * ablwv2 * humid[k])
    end

    # Cloudy layers (free troposphere)
    acloud = cloudc * ablcl2

    for k in n_stratosphere_levels+1:nl1
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
    stratc[1] = stratz * norm_pres
    for k in 2:n_stratosphere_levels
        eps1 = epslw / (σ_levels_thick[k-1] + σ_levels_thick[k])
        stratc[k] = eps1 * norm_pres
    end
end
