"""
Calculate the Legendre weights. TODO reference
"""
function get_legendre_weights(geometry::Geometry)

    @unpack nlat = geometry
    nlat_half = nlat ÷ 2

    z1 = 2.0

    # preallocate the weights
    wt = zeros(nlat_half)

    for i in 1:nlat_half
        local pp
        z = cos(π*(i - 0.25)/(nlat + 0.5))

        while abs(z - z1) > eps(Float64)
            p1 = 1.0
            p2 = 0.0

            for j in 1:nlat
                p3, p2 = p2, p1
                p1 = ((2j - 1.0)*z*p2 - (j - 1.0)*p3)/j
            end

            pp = nlat*(z*p1 - p2)/(z^2 - 1.0)
            z1 = z
            z = z1 - p1/pp
        end

        wt[i] = 2.0/((1.0 - z^2)*pp^2)
    end
    return wt
end

"""
Calculate the Legendre polynomials. TODO reference
"""
function get_legendre_polynomials(  j::Int,
                                    ε::AbstractMatrix,
                                    ε⁻¹::AbstractMatrix,
                                    geometry::Geometry)

    @unpack mx, nx, coslat_half, sinlat_half = geometry

    small = 1e-30
    alp = zeros(mx+1,nx)

    y = coslat_half[j]
    x = sinlat_half[j]

    # Start recursion with n = 1 (m = l) diagonal
    alp[1,1] = sqrt(0.45)
    for m in 2:mx+1
        alp[m,1] = sqrt(0.5*(2.0*(m - 1) + 1.0)/(m - 1))*y*alp[m-1,1]
    end

    # Continue with other elements
    for m in 1:mx+1
        alp[m,2] = (x*alp[m,1])*ε⁻¹[m,2]
    end

    for n in 3:nx
        for m in 1:mx+1
            alp[m,n] = (x*alp[m,n-1] - ε[m,n-1]*alp[m,n-2])*ε⁻¹[m,n]
        end
    end

    # Zero polynomials with absolute values smaller than 10^(-30)
    for n in 1:nx
        for m in 1:mx+1
            if abs(alp[m,n]) <= small
                alp[m,n] = 0.0
            end
        end
    end

    # pick off the required polynomials
    return alp[1:mx,1:nx]
end

"""
Computes the inverse Legendre transform.
"""
function legendre_inverse(  input::Array{T,2},
                            spectral_trans::SpectralTrans{T}
                            geometry::Geometry{T}) where {T<:AbstractFloat}

    @unpack leg_weight, nsh2, leg_poly = spectral_trans
    @unpack nlat, trunc, mx, nx = geometry
    nlat_half = nlat ÷ 2

    # Initialize output array
    output = zeros(Complex{T}, mx, nlat)

    # Loop over Northern Hemisphere, computing odd and even decomposition of incoming field
    for j in 1:nlat_half
        j1 = nlat + 1 - j

        # Initialise arrays TODO preallocate them in a separate struct
        even = zeros(Complex{T}, mx)
        odd  = zeros(Complex{T}, mx)

        # Compute even decomposition
        for n in 1:2:nx
            for m in 1:nsh2[n]
                even[m] = even[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute odd decomposition
        for n in 2:2:nx
            for m in 1:nsh2[n]
                odd[m] = odd[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute Southern Hemisphere
        output[:,j1] = even + odd

        # Compute Northern Hemisphere
        output[:,j]  = even - odd
    end
    return output
end

"""
Computes the Legendre transform
"""
function legendre(  input::Array{T,2},
                    spectral_trans::SpectralTrans{T},
                    geometry::Geometry{T}) where {T<:AbstractFloat}

    @unpack leg_weight, nsh2, leg_poly = spectral_trans
    @unpack nlat, trunc, mx, nx = geometry
    nlat_half = nlat ÷ 2

    # Initialise output array
    output = zeros(Complex{T}, mx, nx)

    even = zeros(Complex{T}, mx, nlat_half)
    odd  = zeros(Complex{T}, mx, nlat_half)

    # Loop over Northern Hemisphere, computing odd and even decomposition of
    # incoming field. The Legendre weights (leg_weight) are applied here
    for j in 1:nlat_half
        # Corresponding Southern Hemisphere latitude
        j1 = nlat + 1 - j

        even[:,j] = (input[:,j1] + input[:,j])*leg_weight[j]
        odd[:,j]  = (input[:,j1] - input[:,j])*leg_weight[j]
    end

    # The parity of an associated Legendre polynomial is the same
    # as the parity of n' - m'. n', m' are the actual total wavenumber and zonal
    # wavenumber, n and m are the indices used for SPEEDY's spectral packing.
    # m' = m - 1 and n' = m + n - 2, therefore n' - m' = n - 1

    # Loop over coefficients corresponding to even associated Legendre polynomials
    for n in 1:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], even[m,:])
        end
    end

    # Loop over coefficients corresponding to odd associated Legendre polynomials
    for n in 2:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], odd[m,:])
        end
    end
    
    return output
end
