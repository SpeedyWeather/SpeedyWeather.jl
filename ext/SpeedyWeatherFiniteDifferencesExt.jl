module SpeedyWeatherFiniteDifferencesExt 

using SpeedyWeather 
import FiniteDifferences
import FiniteDifferences: to_vec 

# FiniteDifferences needs to be able to convert data structures to Vectors and back 
# This doesn't work out of the box with our data types, so we'll define those 
# conversions here.
function FiniteDifferences.to_vec(x::Grid) where Grid <: AbstractGridArray
    x_vec, from_vec = FiniteDifferences.to_vec(Array(x))

    function GridArray_from_vec(x_vec)
        return Grid(reshape(from_vec(x_vec), size(x)), x.nlat_half)
    end 

    return x_vec, GridArray_from_vec
end 

function FiniteDifferences.to_vec(x::LTA) where LTA <: LowerTriangularArray
    x_vec, from_vec = FiniteDifferences.to_vec(x.data)

    function LowerTriangularArray_from_vec(x_vec)
        return LowerTriangularArray(reshape(from_vec(x_vec), size(x)), x.m, x.n)
    end 

    return x_vec, LowerTriangularArray_from_vec
end

end 