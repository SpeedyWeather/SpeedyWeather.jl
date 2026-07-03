module ArrayDimensions

abstract type AbstractArrayDimensions end
abstract type TwoDimensions <: AbstractArrayDimensions end
abstract type ThreeDimensions <: AbstractArrayDimensions end
abstract type FourDimensions <: AbstractArrayDimensions end

# Fields
"""Field Dimensions: 2D longitude-latitude"""
struct XY <: TwoDimensions end
Base.ndims(::XY) = 1

"""Field Dimensions: 3D horizontal + vertical"""
struct XYZ <: ThreeDimensions end
Base.ndims(::XYZ) = 2

"""Field Dimensions: 3D horizontal + time"""
struct XYT <: ThreeDimensions end
Base.ndims(::XYT) = 2

"""4Field Dimensions: D horizontal + vertical + time"""
struct XYZT <: FourDimensions end
Base.ndims(::XYZT) = 3

# ColumnFields
"""ColumnField Dimensions: 3D vertical + horizontal"""
struct ZXY <: ThreeDimensions end
Base.ndims(::ZXY) = 2

"""ColumnField Dimensions: 4D vertical + horizontal + time"""
struct ZXYT <: FourDimensions end
Base.ndims(::ZXYT) = 3

# Spectral
"""Dimension of LowerTriangularArray: 2D degree and order"""
struct LM <: TwoDimensions end
Base.ndims(::LM) = 1

"""Dimension of LowerTriangularArray: 3D horizontal + vertical"""
struct LMZ <: ThreeDimensions end
Base.ndims(::LMZ) = 2

"""Dimension of LowerTriangularArray: 3D horizontal + time"""
struct LMT <: ThreeDimensions end
Base.ndims(::LMT) = 2

"""Dimension of LowerTriangularArray: 4D horizontal + vertical + time"""
struct LMZT <: FourDimensions end
Base.ndims(::LMZT) = 3

default_field_dimensions(::AbstractArray) = XY()                # other dimensions are unspecified
default_column_field_dimensions(::AbstractVector) = XY()
default_column_field_dimensions(::AbstractArray) = ZXY()        # because it's in the name
default_lta_dimensions(::AbstractArray) = LM()                  # other dimensions unspecified

@inline hasvertical(dims) = hasvertical(typeof(dims))
@inline hasvertical(::Type{<:TwoDimensions}) = false
@inline hasvertical(::Type{<:ThreeDimensions}) = true
@inline hasvertical(::Type{XYT}) = false
@inline hasvertical(::Type{LMT}) = false
@inline hasvertical(::Type{<:FourDimensions}) = true

@inline hastime(dims) = hastime(typeof(dims))
@inline hastime(::Type{<:TwoDimensions}) = false
@inline hastime(::Type{<:ThreeDimensions}) = true
@inline hastime(::Type{XYZ}) = false
@inline hastime(::Type{LMZ}) = false
@inline hastime(::Type{<:FourDimensions}) = true

Base.dims2string(dims::AbstractArrayDimensions) = string(nameof(typeof(dims)))

end
