# super type to be used for ZeroBased (0-based), OneBased (1-based) indexing of the spherical harmonics
abstract type IndexBasis end

"""Abstract type to dispatch for 0-based indexing of the spherical harmonic
degree l and order m, i.e. l=m=0 is the mean, the zonal modes are m=0 etc.
This indexing is more common in mathematics."""
abstract type ZeroBased <: IndexBasis end

"""Abstract type to dispatch for 1-based indexing of the spherical harmonic
degree l and order m, i.e. l=m=1 is the mean, the zonal modes are m=1 etc.
This indexing matches Julia's 1-based indexing for arrays."""
abstract type OneBased <: IndexBasis end
