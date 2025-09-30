module LowerTriangularArraysFiniteDifferencesExt 

    import FiniteDifferences: FiniteDifferences, to_vec 
    using LowerTriangularArrays

    function FiniteDifferences.to_vec(x::LTA) where LTA <: LowerTriangularArray
        x_vec, from_vec = FiniteDifferences.to_vec(x.data)

        function LowerTriangularArray_from_vec(x_vec)
            return LowerTriangularArray(reshape(from_vec(x_vec), size(x)), x.spectrum)
        end 

        return x_vec, LowerTriangularArray_from_vec
    end

end 