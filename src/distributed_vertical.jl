for func_name in (:gradient_longitude!, :gradient_latitude!)
    @eval begin
        function $func_name(Out::AbstractArray{NF1,3},
                            In::AbstractArray{NF2,3},
                            args...) where {NF1,NF2}
            
            for k in 1:size(Out)[end]
                Out_layer = view(Out,:,:,k)
                In_layer = view(In,:,:,k)
                $func_name(Out_layer,In_layer,args...)
            end
        end
    end
end