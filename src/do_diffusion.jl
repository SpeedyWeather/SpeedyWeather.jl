"""
    horizontal_diffusion!(  tendency::AbstractMatrix{Complex{NF}}, # tendency of a 
                            A::AbstractMatrix{Complex{NF}},        # spectral horizontal field
                            damp_expl::AbstractMatrix{NF},         # explicit spectral damping
                            damp_impl::AbstractMatrix{NF}          # implicit spectral damping
                            ) where {NF<:AbstractFloat}

Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `damp_expl` and an implicit part `damp_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( tendency::LowerTriangularMatrix{Complex{NF}},   # tendency of a 
                                A::LowerTriangularMatrix{Complex{NF}},          # spectral horizontal field
                                damp_expl::LowerTriangularMatrix{NF},           # explicit spectral damping
                                damp_impl::LowerTriangularMatrix{NF}            # implicit spectral damping
                                ) where {NF<:AbstractFloat}
    
    for lm in eachharmonic(tendency,A,damp_expl,damp_impl)
        @inbounds tendency[lm] = (tendency[lm] - damp_expl[lm]*A[lm])*damp_impl[lm]
    end
end

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::BarotropicModel,
                                lf::Int=1,                          # leapfrog index used
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor = progn.leapfrog[lf]
    @unpack vor_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend,vor,damping,damping_impl)    # diffusion of vorticity
end

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::ShallowWaterModel,
                                lf::Int=1,                          # leapfrog index used
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor, div = progn.leapfrog[lf]
    @unpack vor_tend,div_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend,vor,damping,damping_impl)        # diffusion of vorticity
    horizontal_diffusion!(div_tend,div,damping,damping_impl)        # diffusion of divergence
end

function horizontal_diffusion!( progn::PrognosticVariablesLeapfrog,
                                diagn::DiagnosticVariablesLayer,
                                M::PrimitiveEquationModel,
                                lf::Int=1,                          # leapfrog index used
                                )

    @unpack damping, damping_impl = M.horizontal_diffusion
    @unpack vor,div,temp = progn.leapfrog[lf]
    @unpack vor_tend,div_tend,temp_tend = diagn.tendencies

    horizontal_diffusion!(vor_tend, vor, damping,damping_impl)      # diffusion of vorticity
    horizontal_diffusion!(div_tend, div, damping,damping_impl)      # diffusion of divergence
    horizontal_diffusion!(temp_tend,temp,damping,damping_impl)      # diffusion of temperature
end