# ONE ARGUMENT FUNCTIONS
for func_name in (:scale_coslat!,:scale_coslat²!,:scale_coslat⁻¹!,:scale_coslat⁻²!,
                    :spectral_truncation!,:flipsign!)
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

# TWO ARGUMENT FUNCTIONS
for func_name in (:gradient_longitude!, :gradient_latitude!,
                    :gridded!, :spectral!, :∇⁻²!, :∇²!, :add_tendencies!)
    @eval begin
        function $func_name(Out::AbstractArray{NF1,3},
                            In::AbstractArray{NF2,3},
                            args...;
                            kwargs...) where {NF1,NF2}
            
            for k in 1:size(Out)[end]
                Out_layer = view(Out,:,:,k)
                In_layer = view(In,:,:,k)
                $func_name(Out_layer,In_layer,args...;kwargs...)
            end
        end
    end
end

# THREE ARGUMENT FUNCTIONS
for func_name in (:add_tendencies!,:divergence!,:_divergence!,:curl!)
    @eval begin
        function $func_name(Out::AbstractArray{NF1,3},
                            In1::AbstractArray{NF2,3},
                            In2::AbstractArray{NF2,3},
                            args...;
                            kwargs...) where {NF1,NF2}
            
            for k in 1:size(Out)[end]
                Out_layer = view(Out,:,:,k)
                In1_layer = view(In1,:,:,k)
                In2_layer = view(In2,:,:,k)
                $func_name(Out_layer,In1_layer,In2_layer,args...;kwargs...)
            end
        end
    end
end

# THREE ARGUMENT FUNCTIONS
for func_name in (:bernoulli_potential!,)
    @eval begin
        function $func_name(Out::AbstractArray{NF,3},
                            In1::AbstractArray{NF,3},
                            In2::AbstractArray{NF,3},
                            In3::AbstractArray{NF,2},
                            args...) where NF
            
            for k in 1:size(Out)[end]
                Out_layer = view(Out,:,:,k)
                In1_layer = view(In1,:,:,k)
                In2_layer = view(In2,:,:,k)
                In3_layer = view(In3,:,:)
                $func_name(Out_layer,In1_layer,In2_layer,In3_layer,args...)
            end
        end
    end
end

# FOUR ARGUMENT FUNCTIONS
for func_name in (:UV_from_vordiv!,)
    @eval begin
        function $func_name(Out1::AbstractArray{NF,3},
                            Out2::AbstractArray{NF,3},
                            In1::AbstractArray{NF,3},
                            In2::AbstractArray{NF,3},
                            args...) where NF
            
            for k in 1:size(Out1)[end]
                Out1_layer = view(Out1,:,:,k)
                Out2_layer = view(Out2,:,:,k)
                In1_layer = view(In1,:,:,k)
                In2_layer = view(In2,:,:)
                $func_name(Out1_layer,Out2_layer,In1_layer,In2_layer,args...)
            end
        end
    end
end

# FIVE ARGUMENT FUNCTIONS
for func_name in (:vorticity_fluxes!,)
    @eval begin
        function $func_name(Out1::AbstractArray{NF1,3},
                            Out2::AbstractArray{NF1,3},
                            In1::AbstractArray{NF2,3},
                            In2::AbstractArray{NF2,3},
                            In3::AbstractArray{NF2,3},
                            args...) where {NF1,NF2}
            
            for k in 1:size(Out1)[end]
                Out1_layer = view(Out1,:,:,k)
                Out2_layer = view(Out2,:,:,k)
                In1_layer = view(In1,:,:,k)
                In2_layer = view(In2,:,:,k)
                In3_layer = view(In3,:,:,k)
                $func_name( Out1_layer,
                            Out2_layer,
                            In1_layer,
                            In2_layer,
                            In3_layer,
                            args...)
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