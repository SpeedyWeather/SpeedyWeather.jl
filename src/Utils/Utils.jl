module Utils 

    using DocStringExtensions

    export isincreasing, isdecreasing, clip_negatives!, underflow!
    export flipsign!, nans, print_fields
    export Second, Minute, Hour, Day, Week, second 

    export configure_kernel, launch!

    include("utility_functions.jl")
    include("kernel_launching.jl")
end 