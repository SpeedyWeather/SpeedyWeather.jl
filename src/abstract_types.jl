# MODELS
abstract type ModelSetup end
abstract type Barotropic <: ModelSetup end
abstract type ShallowWater <: ModelSetup end
abstract type PrimitiveEquation <: ModelSetup end
abstract type PrimitiveDry <: PrimitiveEquation end
abstract type PrimitiveWet <: PrimitiveEquation end

abstract type AbstractPlanet end
abstract type AbstractAtmosphere end

# GEOMETRY, GRID
abstract type AbstractGeometry{NF} end
abstract type VerticalCoordinates end

# CONSTANTS (face the dynamical core and not the user)
abstract type AbstractDynamicsConstants{NF} end

# INITIAL CONDITIONS AND OROGRAPHY/BOUNDARIES
abstract type InitialConditions end
abstract type AbstractOrography end

# ATMOSPHERIC COLUMN FOR PARAMETERIZATIONS
abstract type AbstractColumnVariables{NF} end

# FORCING (Barotropic and ShallowWaterModel)
abstract type AbstractForcing{NF} end

# PARAMETERIZATIONS
abstract type AbstractParameterization{NF} end
abstract type BoundaryLayer{NF} <: AbstractParameterization{NF} end
abstract type TemperatureRelaxation{NF} <: AbstractParameterization{NF} end
abstract type VerticalDiffusion{NF} <: AbstractParameterization{NF} end

# INPUT/OUTPUT
abstract type AbstractFeedback end
abstract type AbstractOutputWriter end

# NUMERICS
abstract type HorizontalDiffusion{NF} end
abstract type AbstractImplicit{NF} end
abstract type TimeStepper{NF} end