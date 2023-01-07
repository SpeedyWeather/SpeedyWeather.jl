"""
Check global mean temperature at every vertical level from the m=n=1 coefficient
of the associated Legendre polynomials.
"""
function check_global_mean_temperature(Tabs::Array{Complex{T}, 3},
                                       P::Parameters) where {T}
    @unpack print_dialog = P

    # The m=n=1 coefficent in spectral space is the mean value
    # √2 is due to the normalization of the associated Legendre polynomials
    gmt = sqrt(0.5) * Float64.(real.(Tabs[1, 1, :]))

    if print_dialog
        println(gmt)
    end
end

# # Prints global means of eddy kinetic energy and temperature.
# # Also stops the integration if the computed diagnostics are outside of
# # allowable ranges.
# function check_diagnostics(geometry::Geometry, spectral_trans::SpectralTrans, ξ, D, Tₐ, step,
#                            print_diag=true)
#     @unpack nlev, mx, nx = geometry
#
#     diagnostics = zeros(Float64, nlev, 3)
#
#     # 1. Get global-mean temperature and compute eddy kinetic energy
#     for k in 1:nlev
#         diagnostics[k,3] = √(0.5)*real(Tₐ[1,1,k])
#
#         temp = ∇⁻²(spectral_trans, ξ[:,:,k])
#
#         for m in 2:mx
#             for n in 1:nx
#                 diagnostics[k,1] = diagnostics[k,1] - real(temp[m,n]*conj(ξ[m,n,k]))
#             end
#         end
#
#         temp = ∇⁻²(spectral_trans, D[:,:,k])
#
#         for m in 2:mx
#             for n in 1:nx
#                 diagnostics[k,2] = diagnostics[k,2] - real(temp[m,n]*conj(D[m,n,k]))
#             end
#         end
#     end
#
#     # 2. Print results to screen
#     if print_diag
#         println("Step = $step")
#         println(diagnostics[:,1])
#         println(diagnostics[:,2])
#         println(diagnostics[:,3])
#     end
#
#     # 3. Stop integration if model variables are out of range
#     for k in 1:nlev
#         if diagnostics[k,1] > 500.0 || diagnostics[k,2] > 500.0 ||
#             diagnostics[k,3] < 180.0 || diagnostics[k,3] > 320.0
#
#             println(diagnostics[:,1])
#             println(diagnostics[:,2])
#             println(diagnostics[:,3])
#
#             throw("Model variables out of accepted range")
#         end
#     end
# end
