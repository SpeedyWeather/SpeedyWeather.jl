"""
    function initialise_longwave_radiation!(
        P::Parameters
    )

Initialise variables and parameters used by the longwave radiation parametrization 
"""
function initialise_longwave_radiation!(K::ParameterizationConstants,
                                        P::Parameters) 
    radset!(K,P)
end

"""
    function longwave_radiation!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    )

Compute air temperature tendencies from longwave radiation for an atmospheric column.
For more details see http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf
"""
function longwave_radiation!(   column::ColumnVariables,
                                model::PrimitiveEquation)

    radlw_down!(column, model)
    compute_bbe!(column, model)
    radlw_up!(column, model)
end

"""
    function radset!(model::PrimitiveEquation) where {NF<:AbstractFloat}
        
Compute energy fractions in four longwave bands as a function of temperature.

"""
function radset!(   K::ParameterizationConstants,
                    P::Parameters)

    @unpack NF, nband = P
    @unpack epslw = P.radiation_coefs
    @unpack fband = K

    @assert nband == 4 "Only four bands are supported, given $nband"

    eps1 = 1 - epslw
    for jtemp = 200:320
        fband[jtemp, 2] = (NF(0.148) - NF(3.0e-6) * (jtemp - 247)^2) * eps1
        fband[jtemp, 3] = (NF(0.356) - NF(5.2e-6) * (jtemp - 282)^2) * eps1
        fband[jtemp, 4] = (NF(0.314) + NF(1.0e-5) * (jtemp - 315)^2) * eps1
        fband[jtemp, 1] = eps1 - (fband[jtemp, 2] + fband[jtemp, 3] + fband[jtemp, 4])
    end
    for jb = 1:4
        for jtemp = 100:199
            fband[jtemp, jb] = fband[200, jb]
        end
        for jtemp = 321:400
            fband[jtemp, jb] = fband[320, jb]
        end
    end
end


"""
    function radlw_down!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    ) where {NF<:AbstractFloat}

Compute the downward flux of longwave radiation. Inputs variables are `temp``, `wvi`, `tau2`.
Output column varables are `fsfcd`, `dfabs`, `flux`, `st4a`.
"""
function radlw_down!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}

    @unpack nlev, temp, wvi, tau2 = column
    @unpack nband, sbc = model.parameters
    @unpack n_stratosphere_levels = model.geometry
    @unpack fband = model.parameterization_constants
    @unpack epslw, emisfc = model.parameters.radiation_coefs
    
    # 1. Blackbody emission from atmospheric levels.
    #    The linearized gradient of the blakbody emission is computed
    #    from temperatures at layer boundaries, which are interpolated 
    #    assuming a linear dependence of T on log_sigma.
    #    Above the first (top) level, the atmosphere is assumed isothermal.
    # 
    # Temperature at level boundaries
    for k=1:nlev-1
        column.st4a[k, 1] = temp[k] + wvi[k, 2] * (temp[k + 1] - temp[k])
    end
    
    # Mean temperature in stratospheric layers
    for k=1:n_stratosphere_levels-1
        column.st4a[k, 2] = NF(0.75) * temp[k] + NF(0.25) * column.st4a[k, 1]
        column.st4a[k+1, 2] = NF(0.50) * temp[k+1] + NF(0.25) * (column.st4a[k, 1] + column.st4a[k+1, 1])
    end

    # Temperature gradient in tropospheric layers
    anis = NF(1.)
    anish = NF(0.5) # 0.5 * anis

    for k = n_stratosphere_levels+1:nlev-1
        column.st4a[k, 2] = anish * max(column.st4a[k, 1] - column.st4a[k - 1, 1], NF(0.))
    end
    column.st4a[nlev, 2] = anis * max(temp[nlev] - column.st4a[nlev-1, 1], NF(0.))

    # Blackbody emission in the stratosphere
    for k = 1:n_stratosphere_levels
        column.st4a[k, 1] = sbc * column.st4a[k, 2]^4
        column.st4a[k, 2] = NF(0.)
    end

    # Blackbody emission in the troposphere
    for k = n_stratosphere_levels+1:nlev
            st3a = sbc * temp[k]^3
            column.st4a[k, 1] = st3a * temp[k]
            column.st4a[k, 2] = NF(4.) * st3a * column.st4a[k, 2]
    end

    # 2.0 Initialization of fluxes
    column.fsfcd = 0.
    column.dfabs .= 0.

    # 3.0 Emission ad absorption of longwave downward flux.
    #    For downward emission, a correction term depending on the
    #    local temperature gradient and on the layer transmissivity is
    #    added to the average (full - level) emission of each layer.

    # 3.1  Stratosphere
    k = n_stratosphere_levels-1
    for jb = 1:nband-2
        emis = 1. - tau2[k, jb]
        brad = fband[convert(Int, temp[k]), jb] * (column.st4a[k, 1] + emis * column.st4a[k, 2])
        column.flux[jb] = emis * brad
        column.dfabs[k] = column.dfabs[k] - column.flux[jb]
    end
    column.flux[nband-1:nband] .= NF(0.)

    # 3.2  Troposphere
    for jb = 1:nband
        for k = n_stratosphere_levels:nlev
            emis = 1. - tau2[k, jb]
            brad = fband[convert(Int, temp[k]), jb] * (column.st4a[k, 1] + emis * column.st4a[k, 2])
            column.dfabs[k] = column.dfabs[k] + column.flux[jb]
            column.flux[jb] = tau2[k, jb] * column.flux[jb] + emis * brad
            column.dfabs[k] = column.dfabs[k] - column.flux[jb]
        end
    end

    # 3.3 Surface downward flux
    for jb = 1: nband
        column.fsfcd = column.fsfcd + emisfc * column.flux[jb]
    end
    eps1 = epslw * emisfc
    
    # 3.4 Correction for "black" band (incl. surface reflection)
    corlw = eps1 * column.st4a[nlev, 1]
    column.dfabs[nlev] = column.dfabs[nlev] - corlw
    column.fsfcd = column.fsfcd + corlw
end

"""
    function compute_bbe!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    ) where {NF<:AbstractFloat}

# Computes black-body (or grey-body) emissions.
Input and output variables are `ts` and `fsfcu`, respectively.
"""
function compute_bbe!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}

@unpack sbc = model.parameters
@unpack emisfc = model.parameters.radiation_coefs
@unpack ts = column

    column.fsfcu = emisfc * sbc * ts^4
end

"""
    function radlw_up!(
        column::ColumnVariables{NF}, model::PrimitiveEquation
    ) where {NF<:AbstractFloat}

# Computes the upward flux of longwave radiation.
Input variables are  `nlev`, `temp`, `fsfcu`, `fsfcd`, `flux`, `ts`, `tau2`, `st4a`, `dfabs`, `stratc`, `σ_levels_thick`, `n_stratosphere_levels`. Output column variables are `fsfc` and `ftop`.
"""
function radlw_up!(
    column::ColumnVariables{NF}, model::PrimitiveEquation
) where {NF<:AbstractFloat}

@unpack nband = model.parameters
@unpack epslw, emisfc = model.parameters.radiation_coefs
@unpack fband = model.parameterization_constants
@unpack nlev, temp, fsfcu, fsfcd, flux, ts, tau2, st4a, dfabs, stratc = column
@unpack σ_levels_thick, n_stratosphere_levels = model.geometry

    column.fsfc = fsfcu - fsfcd
    
    refsfc = 1 - emisfc
    for jb = 1: nband
        column.flux[jb] = fband[convert(Int, ts), jb] * fsfcu + refsfc * flux[jb]
    end

    # 4.2  Troposphere
    # Correction for "black" band
    dfabs[nlev] = dfabs[nlev] + epslw * fsfcu
    for jb = 1: nband
        for k = nlev:-1:n_stratosphere_levels
            emis = 1. - tau2[k, jb]
            brad = fband[convert(Int, temp[k]), jb] * (st4a[k, 1] - emis * st4a[k, 2])
            column.dfabs[k] = dfabs[k] + flux[jb]
            column.flux[jb] = tau2[k, jb] * flux[jb] + emis * brad
            column.dfabs[k] = dfabs[k] - flux[jb]
        end
    end
    # 4.3  Stratosphere
    k = n_stratosphere_levels - 1
    for jb = 1: 2
        emis = 1 - tau2[k, jb]
        brad = fband[convert(Int, temp[k]), jb] * (st4a[k, 1] - emis * st4a[k, 2])
        column.dfabs[k] = dfabs[k] + flux[jb]
        column.flux[jb] = tau2[k, jb] * flux[jb] + emis * brad
        column.dfabs[k] = dfabs[k] - flux[jb]
    end

    # Correction for "black" band and polar night cooling
    column.ftop = 0. # Init
    for k=1:n_stratosphere_levels-1
        corlw1 = σ_levels_thick[k] * stratc[k+1] * st4a[k, 1] + stratc[k]
        corlw2 = σ_levels_thick[k+1] * stratc[k+1] * st4a[k+1, 1]
        dfabs[k] -= corlw1
        dfabs[k+1] -= corlw2
        column.ftop += corlw1 + corlw2
    end

    # 4.4  Outgoing longwave radiation
    for jb = 1:nband
        column.ftop = column.ftop + flux[jb]
    end
end
