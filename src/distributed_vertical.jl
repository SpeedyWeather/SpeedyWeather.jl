# TWO ARGUMENT FUNCTIONS
for func_name in (:gradient_longitude!, :gradient_latitude!,
                    :gridded!, :spectral!, :∇⁻²!)
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

# ONE ARGUMENT FUNCTIONS
for func_name in (:unscale_coslat!,:scale_coslat!,
                    :spectral_truncation!)
    @eval begin
        function $func_name(A::AbstractArray{NF,3},
                            args...) where NF

            for k in 1:size(A)[end]
                A_layer = view(A,:,:,k)
                $func_name(A_layer,args...)
            end
        end
    end
end

"""Leapfrog! for 3D arrays that loops over all vertical layers."""
function leapfrog!( A::AbstractArray{Complex{NF},4},        # a prognostic variable (spectral)
                    tendency::AbstractArray{Complex{NF},3}, # tendency (dynamics + physics) of A
                    dt::Real,                               # time step (mostly =2Δt, but for init steps =Δt,Δt/2)
                    C::Constants{NF},                       # struct containing all constants used at runtime
                    lf::Int=2                               # leapfrog index to dis/enable(default) William's filter
                    ) where {NF<:AbstractFloat}             # number format NF

    for k in 1:size(A)[end]                         # loop over vertical levels (last dimension)
        A_layer = view(A,:,:,:,k)                   # extract vertical layers as views to not allocate any memory
        tendency_layer = view(tendency,:,:,k)
        leapfrog!(A_layer,tendency_layer,dt,C,lf)   # make a timestep forward for each layer
    end
end