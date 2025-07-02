module Utils 

using DocStringExtensions

    import Dates: Dates, DateTime, Period, Millisecond, Second, Minute, Hour, Day, Week, Month, Year

    export isincreasing, isdecreasing, clip_negatives!, underflow!
    export flipsign!, nans, print_fields

    export DateTime, Millisecond, Second, Minute, Hour, Day, Week, Month, Year, Century, Millenium

    export configure_kernel, launch!

    include("utility_functions.jl")
    include("kernel_launching.jl")
end 