module SpeedyWeatherInternals

    # ARCHITECTURES (Device handling)
    include("Architectures/Architectures.jl")

    # UTILS (Kernel launching and various utilities)
    include("Utils/Utils.jl")

    # SPEEDY PARAMETERS (parameter handling)
    include("SpeedyParameters/SpeedyParameters.jl")

end
