struct PhyPar{T<:AbstractFloat}
    # SIZE OF SPECTRAL GRID
    trunc::Int      # Spectral truncation
    nx::Int         # Number of total wavenumbers
    mx::Int         # Number of zonal wavenumbers

    # LEGENDRE ARRAYS
    leg_weight::Array{T,1}          # Legendre weights
    nsh2::Array{Int,1}              # What's this?
    leg_poly::Array{Complex{T},3}   # Legendre polynomials

    # HARMONIC AND BIHARMONIC DIFFUSION
    ∇²::Array{T,2}          # Laplacian = l*(l+1)/(R_earth^2)
    ∇⁻²::Array{T,2}         # inverse Laplacian
    ∇⁴::Array{T,2}          # Laplacian squared, for biharmonic diffusion

    # Quantities required by functions grad, uvspec, and vds
    gradx::Array{T,1}
    uvdx::Array{T,2}
    uvdym::Array{T,2}
    uvdyp::Array{T,2}
    gradym::Array{T,2}
    gradyp::Array{T,2}
    vddym::Array{T,2}
    vddyp::Array{T,2}
end

"""
Compute physical parametrization tendencies
"""
function phypar{T}(utend, vtend, ttend, trtend) where {T<:AbstractFloat}


end
