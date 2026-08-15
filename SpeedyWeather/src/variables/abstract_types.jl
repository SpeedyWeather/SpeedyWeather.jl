abstract type AbstractVariableDim end                           # Dimensions of variables

abstract type AbstractVariableDim2D <: AbstractVariableDim end  # 2D group XY, LM
abstract type AbstractVariableDim3D <: AbstractVariableDim end  # 3D group XYZ, XYT, LMZ, LMT
abstract type AbstractVariableDim4D <: AbstractVariableDim end  # 4D group XYZT, LMZT

abstract type AbstractVariable{D <: AbstractVariableDim} end    # Variable type with dimensions D
abstract type AbstractVariables end                             # Container for all variable arrays in the model