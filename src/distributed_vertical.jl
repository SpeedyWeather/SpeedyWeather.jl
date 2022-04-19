for func_name in (:gradient_longitude!, :gradient_latitude!,
                    :gridded!, :spectral!)
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

for func_name in (:unscale_coslat!,:scale_coslat!)
    @eval begin
        function $func_name(A::AbstractArray{NF,3},
                            G::Geometry{NF}) where NF

            for k in 1:size(A)[end]
                A_layer = view(A,:,:,k)
                $func_name(A_layer,G)
            end
        end
    end
end