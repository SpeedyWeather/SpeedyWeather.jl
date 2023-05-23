struct NoForcing{NF} <: AbstractForcing{NF} end
NoForcing(SG::SpectralGrid) = NoForcing{SG.NF}()

function initialize!(   forcing::NoForcing,
                        model::ModelSetup)
    return nothing
end

function forcing!(  diagn::DiagnosticVariablesLayer,
                    forcing::NoForcing,
                    time::DateTime)
    return nothing
end

# function interface_relaxation!( η::LowerTriangularMatrix{Complex{NF}},
#                                 surface::SurfaceVariables{NF},
#                                 time::DateTime,         # time of relaxation
#                                 M::ShallowWater,   # contains η⁰, which η is relaxed to
#                                 ) where NF    

#     (; pres_tend ) = surface
#     (; seasonal_cycle, equinox, axial_tilt ) = M.parameters.planet
#     A = M.parameters.interface_relax_amplitude

#     s = 45/23.5     # heuristic conversion to Legendre polynomials
#     θ = seasonal_cycle ? s*axial_tilt*sin(Dates.days(time - equinox)/365.25*2π) : 0
#     η2 = convert(NF,A*(2sind(θ)))           # l=1,m=0 harmonic
#     η3 = convert(NF,A*(0.2-1.5cosd(θ)))     # l=2,m=0 harmonic

#     τ⁻¹ = inv(M.constants.interface_relax_time)
#     pres_tend[2] += τ⁻¹*(η2-η[2])
#     pres_tend[3] += τ⁻¹*(η3-η[3])
# end