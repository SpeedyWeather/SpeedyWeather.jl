# utils for combining LowerTriangularArrays and Fields for the spectral transforms
# created originally as a workaround for direct allocations for Reactant

"""
$(TYPEDSIGNATURES)
Create a Field similar to the `spec` (similar to underlying `spec.data` field) 
as an output field for spectral transforms.
"""
function Base.similar(spec::LowerTriangularArray, grid::AbstractGrid, T::DataType=eltype(spec.data))
    return Field(similar(spec.data, T, RingGrids.get_npoints(grid), size(spec.data)[2:end]...), grid)
end 

"""
$(TYPEDSIGNATURES)
Create a LowerTriangularArray similar to the `field` (similar to underlying `field.data` field) 
as an output field for spectral transforms.
"""
function Base.similar(field::AbstractField, spectrum::AbstractSpectrum, T::DataType=eltype(field.data))
    return LowerTriangularArray(similar(field.data, T, LowerTriangularArrays.nonzeros(spectrum), size(field.data)[2:end]...), spectrum)
end