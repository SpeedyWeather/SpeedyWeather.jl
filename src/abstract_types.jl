abstract type ModelSetup end
abstract type Barotropic <: ModelSetup end
abstract type ShallowWater <: ModelSetup end
abstract type PrimitiveEquation <: ModelSetup end
abstract type PrimitiveDryCore <: PrimitiveEquation end
abstract type PrimitiveWetCore <: PrimitiveEquation end

abstract type AbstractParameters{M} end

abstract type AbstractOrography{NF} end
abstract type AbstractBoundaries{NF} end

abstract type AbstractColumnVariables{NF} end
abstract type AbstractConstants{NF} end
