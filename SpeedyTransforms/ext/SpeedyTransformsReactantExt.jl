module SpeedyTransformsReactantExt

using SpeedyTransforms
using Reactant

using SpeedyWeatherInternals.Architectures: ReactantDevice

import SpeedyTransforms: AbstractSpectralTransform
import Base: eltype

# TODO::maybe change the reactant logic so that SpectralGrid directly has this NF? but this may cause trouble elsewhere
Base.eltype(::AbstractSpectralTransform{NF, <:ReactantDevice}) where {NF} = Reactant.TracedRNumber{NF}

end
