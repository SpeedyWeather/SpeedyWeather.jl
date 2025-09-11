module RingGridsFiniteDifferencesExt 

    import FiniteDifferences: FiniteDifferences, to_vec 

    function FiniteDifferences.to_vec(x::AbstractField)
        x_vec, from_vec = FiniteDifferences.to_vec(Array(x))
    
        function field_from_vec(x_vec)
            return Field(reshape(from_vec(x_vec), size(x)), x.grid)
        end 
    
        return x_vec, field_from_vec
    end 
    
end 
