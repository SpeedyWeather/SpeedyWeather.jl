module SpeedyParameters

    # STRUCTURE
    using DocStringExtensions

    # ARRAYS
    import ComponentArrays: ComponentArray, ComponentVector, labels, label2index, getaxes

    # UTILITIES
    import MacroTools
    import ModelParameters: ModelParameters, AbstractParam
    import ConstructionBase: constructorof, getproperties, setproperties

    # DOMAINS
    import DomainSets: Domain, RealLine, NonnegativeRealLine, PositiveRealLine, NegativeRealLine, UnitInterval
    using DomainSets.IntervalSets
    export Unbounded, Positive, Nonnegative, Negative
    export ComponentVector

    # Domain aliases
    const Unbounded = RealLine()
    const Positive = PositiveRealLine()
    const Nonnegative = NonnegativeRealLine()
    const Negative = NegativeRealLine()

    # Parameter utilities
    export SpeedyParam, SpeedyParams
    export @parameterized, parameters, parameterof, reconstruct, stripparams, attributes, bounds, description, value
    
    include("parameters.jl")
    
end 