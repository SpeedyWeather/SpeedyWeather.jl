"""
    m = roundup_fft(n::Int)

Returns an integer `m >= n` with only prime factors 2 and 3 to obtain an efficiently
fourier-transformable number of longitudes, m = 2^i * 3^j >= n, with i>=0 but j = 0,1.
"""
function roundup_fft(n::Int)
    lz = leading_zeros(n)                   # determine scale of n

    # create a mask that's 1 for all but the two most significant figures of n
    # for finding the next largest integer of n with factors 2 and 3
    mask = (1 << (8*sizeof(n)-lz-2))-1    
    
    # round up by adding mask as offset and then masking all insignificant bits
    n_roundedup = (n + mask) & ~mask
    return n_roundedup
end

"""
    triangular_truncation(trunc::Int,nlon::Int,nlat::Int)

Tests whether the inputs `trunc, nlon, nlat` satisfy the triangular truncation constraints.
`trunc` is the maximum degree and order of the Legendre polynomials (0-based), `nlon` is the
number of longitudes, `nlat` the number of latitudes on the spatial grid. The constraints are

    - nlon >= 3T+1
    - nlat >= (3T+1)/2
    - nlon = 2nlat."""
function triangular_truncation(trunc::Int,nlon::Int,nlat::Int)
    info = "$(nlon)x$(nlat) grid and T$trunc violate triangular truncation constraints, "
    constraints = "nlon >= 3T+1, nlat >= (3T+1)/2, nlon = 2nlat."
    feedback = info*constraints
    @assert nlon >= 3*trunc+1 feedback
    @assert nlat >= (3*trunc+1)/2 feedback
    @assert nlon == 2nlat feedback
end

"""
    nlon, nlat = triangular_truncation(trunc::Int)

Returns the grid size `nlon,nlat` for a spectral truncation `trunc` that satisfies the
triangular truncation constraints. Returns the smallest pair `nlon,nlat` that is also
easily Fast Fourier-transformable, as determined in `roundup_fft`."""
function triangular_truncation(trunc::Int)
    nlon = roundup_fft(3*trunc+1)
    nlat = nlon÷2
    triangular_truncation(trunc,nlon,nlat)
    return nlon,nlat
end

"""
    trunc = triangular_truncation(nlon::Int,nlat::Int)

Returns the largest spectral truncation `trunc` that satisfies the triangular truncation constraints
based on the grid size `nlon,nlat`, which may or may not be easily Fast Fourier-transformable."""
function triangular_truncation(nlon::Int,nlat::Int)
    trunc = floor(Int,(nlon-1)/3)
    triangular_truncation(trunc,nlon,nlat)
    return trunc
end

"""
    spectral_truncation!(alms,trunc)

Truncate spectral coefficients `alms` in-place by setting (a) the upper right triangle to zero and (b)
all coefficients for which the degree l is larger than the truncation `trunc`."""
function spectral_truncation!(  alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                trunc::Int                          # truncate to total wave number `trunc`
                                ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m of the legendre polynomials

    @inbounds for m in 1:mmax+1                 # order m = 0,mmax but 1-based
        for l in 1:lmax+1                       # degree l = 0,lmax but 1-based
            if m > l || l > trunc+1             # zero upper triangle (m>l) and degrees l>trunc
                alms[l,m] = zero(Complex{NF})
            end
        end
    end
    return alms
end

"""
    spectral_truncation!(alms)

Truncate spectral coefficients `alms` in-place by setting the upper right triangle to zero. This is
to enforce that all coefficients for which the degree l is larger than order m are zero."""
spectral_truncation!(alms::AbstractMatrix) = spectral_truncation!(alms,size(alms)[1])

"""
    alms_trunc = spectral_truncation(alms,trunc)

Returns a spectral coefficient matrix `alms_trunc` that is truncated from `alms` to the size (`trunc`+1)^2.
`alms_trunc` only contains those coefficient of `alms` for which m,l <= trunc, and m>l are zero anyway. If
`trunc` is larger than the implicit truncation in `alms` obtained from its size than `spectral_interpolation`
is automatically called instead, returning `alms_interp`, a coefficient matrix that is larger than `alms`
with padded zero coefficients."""
function spectral_truncation(   alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                trunc::Int                          # truncate to degree and order trunc
                                ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m or the legendre polynomials
    @boundscheck lmax == mmax || throw(BoundsError)
    trunc > lmax && return spectral_interpolation(alms,trunc)

    alms_trunc = Matrix{Complex{NF}}(undef,trunc+1,trunc+1)
    copyto!(alms_trunc,@view(alms[1:trunc+1,1:trunc+1]))
    spectral_truncation!(alms_trunc,trunc)
    return alms_trunc
end

"""
    alms_interp = spectral_interp(alms,trunc)

Returns a spectral coefficient matrix `alms_interp` that is `alms` padded with zeros to interpolate in
spectral space. If `trunc` is smaller or equal to the implicit truncation in `alms` obtained from its size
than `spectral_truncation` is automatically called instead, returning `alms_trunc`, a coefficient matrix that
is smaller than `alms`, implicitly setting higher degrees and orders to zero."""
function spectral_interpolation(    alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                    trunc::Int                          # truncate to degree and order trunc
                                    ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m or the legendre polynomials
    @boundscheck lmax == mmax || throw(BoundsError)
    trunc <= lmax && return spectral_truncation(alms,trunc)

    alms_interp = zeros(Complex{NF},trunc+1,trunc+1)
    copyto!(@view(alms_interp[1:lmax+1,1:mmax+1]),alms)
    return alms_interp
end

# """Spectral truncation with unpacking Geospectral struct."""
# function spectral_truncation!(  A::AbstractArray{Complex{NF},2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral_truncation!(A,G.spectral.trunc)    # unpack GeoSpectral struct
# end

# """Spectral truncation of a grid-point field with memory allocation."""
# function spectral_truncation(   input::AbstractArray{NF,2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     input_spectral = spectral(input,G)          # allocates memory
#     spectral_truncation!(input_spectral,G)      # in-place truncation

#     # allocates memory to return spectrally truncated gridded field
#     return gridded(input_spectral, G)       
# end

# """In-place version of spectral trunction of a grid-point field."""
# function spectral_truncation!(  input::AbstractArray{NF,2},
#                                 input_spectral::AbstractArray{Complex{NF},2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral!(input_spectral,input,G)       # in-place spectral transform from input to input_spectral
#     spectral_truncation!(input_spectral,G)  # in-place truncation
#     gridded!(input,input_spectral,G)        # in-place backtransform
# end

# function spectral_truncation!(  input::Array{NF,2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral_truncation!(input,spectral(input,G),G)
# end