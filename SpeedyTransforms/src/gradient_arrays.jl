"""$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1))."""
recursion_factor(l, m) = sqrt((l^2 - m^2) / (4 * l^2 - 1))

"""$(TYPEDSIGNATURES)     
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` (1-based) of the `spectrum` in number format `NF`."""
function recursion_factors(
        ::Type{NF},             # number format NF
        spectrum::Spectrum,
    ) where {NF}

    ϵ_lm = zeros(NF, spectrum)                              # store in lower triangular matrix
    for (m, degrees) in enumerate(orders(spectrum))         # loop over 1-based l, m
        for l in degrees
            ϵ_lm[l, m] = recursion_factor(l - 1, m - 1)     # convert to 0-based l, m for function arguments
        end
    end
    return ϵ_lm
end

# if number format not provided use DEFAULT_NF
recursion_factors(spectrum::Spectrum) = recursion_factors(DEFAULT_NF, spectrum)

"""$(TYPEDSIGNATURES)
Precomputed factors used for meridional gradients.
On unit sphere, hence 1/radius-scaling is omitted."""
function meridional_gradient_factors(
        ::Type{NF},             # number format NF
        spectrum::Spectrum
    ) where {NF}

    # get recursion factors for spherical harmonics
    (; lmax, mmax) = spectrum
    ϵ_lm = recursion_factors(NF, Spectrum(lmax + 1, mmax))

    # meridional gradient for scalars (coslat scaling included)
    grad_y1 = zeros(NF, spectrum)                           # term 1, mul with harmonic l-1, m
    grad_y2 = zeros(NF, spectrum)                           # term 2, mul with harmonic l+1, m

    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            grad_y1[l, m] = -(l - 2) * ϵ_lm[l, m]
            grad_y2[l, m] = (l + 1) * ϵ_lm[l + 1, m]
        end
        # explicitly set the last row to zero, so that kernels yield correct gradient in last row
        grad_y2[degrees[end], m] = 0
    end

    # meridional gradient used to get from u, v/coslat to vorticity and divergence
    grad_y_vordiv1 = zeros(NF, spectrum)                    # term 1, mul with harmonic l-1, m
    grad_y_vordiv2 = zeros(NF, spectrum)                    # term 2, mul with harmonic l+1, m

    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            grad_y_vordiv1[l, m] = l * ϵ_lm[l, m]
            grad_y_vordiv2[l, m] = (l - 1) * ϵ_lm[l + 1, m]
        end
        # explicitly set the last row to zero, so that kernels yield correct gradient in last row
        grad_y_vordiv1[degrees[end], m] = 0
        grad_y_vordiv2[degrees[end], m] = 0
    end

    return grad_y1, grad_y2, grad_y_vordiv1, grad_y_vordiv2
end

# if number format not provided use default
meridional_gradient_factors(spectrum::Spectrum) = meridional_gradient_factors(DEFAULT_NF, spectrum)

# zonal gradient used to get from u, v/coslat to vorticity and divergence
function zonal_gradient_factors(
        spectrum::Spectrum,
    )
    grad_x_vordiv = zeros(Int, spectrum)                    # allocate as int, can be converted later
    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            grad_x_vordiv[l, m] = m - 1                       # just the 0-based order (=zonal wavenumber), omit imaginary unit
        end # last row zero to get vor and div correct
        grad_x_vordiv[degrees[end], m] = 0
    end

    return grad_x_vordiv
end

# zonal integration (sort of) to get from vorticity and divergence to u, v*coslat
function zonal_integration_factors(
        ::Type{NF},             # number format NF
        spectrum::Spectrum,
    ) where {NF}
    vordiv_to_uv_x = zeros(NF, spectrum)
    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            # -m/(l*(l+1)) is the 0-based integration factor, 1-based with m->m-1, l->l-1
            vordiv_to_uv_x[l, m] = -(m - 1) / (l * (l - 1))
        end
    end
    vordiv_to_uv_x[1] = 0                                   # remove NaN from 0/0
    return vordiv_to_uv_x
end

# if number format not provided use default
zonal_integration_factors(spectrum::Spectrum) = zonal_integration_factors(DEFAULT_NF, spectrum)

# meridional integration (sort of) to get from vorticity and divergence to u, v*coslat
function meridional_integration_factors(
        ::Type{NF},             # number format NF
        spectrum::Spectrum,
    ) where {NF}

    # get recursion factors for spherical harmonics
    (; lmax, mmax) = spectrum
    ϵ_lm = recursion_factors(NF, Spectrum(lmax + 1, mmax))

    vordiv_to_uv1 = zeros(NF, spectrum)                     # term 1, to be mul with harmonic l-1, m
    vordiv_to_uv2 = zeros(NF, spectrum)                     # term 2, to be mul with harmonic l+1, m

    for (m, degrees) in enumerate(orders(spectrum))         # 1-based degree l, order m
        for l in degrees
            vordiv_to_uv1[l, m] = ϵ_lm[l, m] / (l - 1)
            vordiv_to_uv2[l, m] = ϵ_lm[l + 1, m] / l
        end
        # explicitly set the last row of vordiv_to_uv2 to zero, so that kernels yield correct gradient in last row
        vordiv_to_uv2[degrees[end], m] = 0
    end

    vordiv_to_uv1[1] = 0                 # remove NaN from 0/0
    return vordiv_to_uv1, vordiv_to_uv2
end

# if number format not provided use default
meridional_integration_factors(spectrum::Spectrum) = meridional_integration_factors(DEFAULT_NF, spectrum)

"""$(TYPEDSIGNATURES)
Get the eigenvalues of the spherical harmonics up to degree `lmax` and order `mmax`
as determined by the `spectrum` in number format `NF`."""
function get_eigenvalues(::Type{NF}, spectrum::Spectrum) where {NF}
    (; lmax, mmax) = spectrum
    return eigenvalues = [NF(-l * (l + 1)) for l in 0:(lmax - 1)]      # 0:lmax-1 for 0-based indexing
end

get_eigenvalues(spectrum::Spectrum) = get_eigenvalues(DEFAULT_NF, spectrum)

gradient_arrays(spectrum::Spectrum) = gradient_arrays(DEFAULT_NF, spectrum)

"""$(TYPEDSIGNATURES)
Precompute all gradient-related arrays for spherical harmonics defined by `spectrum` in number format `NF`."""
function gradient_arrays(::Type{NF}, spectrum::Spectrum) where {NF}

    grad_y1, grad_y2, grad_y_vordiv1, grad_y_vordiv2 = meridional_gradient_factors(NF, spectrum)
    grad_x_vordiv = zonal_gradient_factors(spectrum)

    # VORTICITY, DIVERGENCE TO U, V includes some precomputed "integration" factors
    vordiv_to_uv1, vordiv_to_uv2 = meridional_integration_factors(NF, spectrum)
    vordiv_to_uv_x = zonal_integration_factors(NF, spectrum)

    eigenvalues = get_eigenvalues(NF, spectrum)     # = -l*(l+1), degree l of spherical harmonic
    eigenvalues⁻¹ = inv.(eigenvalues)
    eigenvalues⁻¹[1] = 0                            # set the integration constant to 0

    gradients = (;
        # GRADIENTS
        grad_y1,
        grad_y2,
        grad_y_vordiv1,
        grad_y_vordiv2,
        grad_x_vordiv,
        # VORTICITY, DIVERGENCE TO U, V
        vordiv_to_uv1,
        vordiv_to_uv2,
        vordiv_to_uv_x,
        # EIGENVALUES
        eigenvalues,
        eigenvalues⁻¹,
    )

    return gradients
end
