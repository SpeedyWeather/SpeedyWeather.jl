module ArrayDimensions

abstract type AbstractArrayDimensions end

"""Super type for all (physical) 2D dimensions, i.e. 1 array dimensions (horizontal is unravelled)"""
abstract type TwoDimensions <: AbstractArrayDimensions end
"""Super type for all (physical) 3D dimensions, i.e. 2 array dimensions (horizontal is unravelled)"""
abstract type ThreeDimensions <: AbstractArrayDimensions end
"""Super type for all (physical) 4D dimensions, i.e. 3 array dimensions (horizontal is unravelled)"""
abstract type FourDimensions <: AbstractArrayDimensions end

# when viewing a getindex dimensions have to be converted, generalize this with a
@inline Base.getindex(dims::AbstractArrayDimensions, args...) = dims

# Fields
"""Field Dimensions: 2D longitude-latitude"""
struct XY <: TwoDimensions end
@inline Base.ndims(::XY) = 1

"""Field Dimensions: 3D horizontal + vertical"""
struct XYZ <: ThreeDimensions end
@inline Base.ndims(::XYZ) = 2
@inline Base.getindex(::XYZ, ::Colon, k::Union{Integer, CartesianIndex}, args...) = XY()

"""Field Dimensions: 3D horizontal + time"""
struct XYT <: ThreeDimensions end
@inline Base.ndims(::XYT) = 2
@inline Base.getindex(::XYT, ::Colon, s::Union{Integer, CartesianIndex}, args...) = XY()

"""Field Dimensions: 4D horizontal + vertical + time"""
struct XYZT <: FourDimensions end
@inline Base.ndims(::XYZT) = 3
@inline Base.getindex(::XYZT, ::Colon, k::Union{Integer, CartesianIndex}, args...) = XYT()
@inline Base.getindex(::XYZT, ::Colon, k, s::Union{Integer, CartesianIndex}, args...) = XYZ()
@inline Base.getindex(::XYZT, ::Colon, k::Union{Integer, CartesianIndex}, s::Union{Integer, CartesianIndex}, args...) = XY() 

# ColumnFields
"""ColumnField Dimensions: 3D vertical + horizontal"""
struct ZXY <: ThreeDimensions end
@inline Base.ndims(::ZXY) = 2

"""ColumnField Dimensions: 4D vertical + horizontal + time"""
struct ZXYT <: FourDimensions end
@inline Base.ndims(::ZXYT) = 3

# Spectral
"""Dimension of LowerTriangularArray: 2D degree and order"""
struct LM <: TwoDimensions end
@inline Base.ndims(::LM) = 1

"""Dimension of LowerTriangularArray: 3D horizontal + vertical"""
struct LMZ <: ThreeDimensions end
@inline Base.ndims(::LMZ) = 2
@inline Base.getindex(::LMZ, ::Colon, k::Union{Integer, CartesianIndex}, args...) = LM()

"""Dimension of LowerTriangularArray: 3D horizontal + time"""
struct LMT <: ThreeDimensions end
@inline Base.ndims(::LMT) = 2
@inline Base.getindex(::LMT, ::Colon, s::Union{Integer, CartesianIndex}, args...) = LM()

"""Dimension of LowerTriangularArray: 4D horizontal + vertical + time"""
struct LMZT <: FourDimensions end
@inline Base.ndims(::LMZT) = 3
@inline Base.getindex(::LMZT, ::Colon, k::Union{Integer, CartesianIndex}, args...) = LMT()
@inline Base.getindex(::LMZT, ::Colon, k, s::Union{Integer, CartesianIndex}, args...) = LMZ()
@inline Base.getindex(::LMZT, ::Colon, k::Union{Integer, CartesianIndex}, s::Union{Integer, CartesianIndex}, args...) = LM()

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

# Dimension grouping unions for dispatch
"""Union of all 2D dimension types"""
const Dimensions2D = Union{XY, LM}
"""Union of all 3D dimension types"""
const Dimensions3D = Union{XYZ, XYT, ZXY, LMZ, LMT}
"""Union of all 4D dimension types"""
const Dimensions4D = Union{XYZT, ZXYT, LMZT}
"""Union of dimension types that have a time dimension"""
const DimensionsWithTime = Union{XYT, XYZT, LMT, LMZT, ZXYT}
"""Union of dimension types that have a vertical dimension"""
const DimensionsWithVertical = Union{XYZ, XYZT, LMZ, LMZT, ZXY, ZXYT}
"""Union of dimension types that have both time and vertical dimensions"""
const DimensionsWithTimeAndVertical = Union{XYZT, ZXYT, LMZT}

Base.dims2string(dims::AbstractArrayDimensions) = string(nameof(typeof(dims)))

# move here as otherwise Method overwritten if
# defined in both LowerTriangularArrays and RingGrids
function Base.DimensionMismatch(data::AbstractArray, dims::AbstractArrayDimensions)
    return DimensionMismatch("Dimensionality of $(summary(data)) does not match $dims")
end

end
