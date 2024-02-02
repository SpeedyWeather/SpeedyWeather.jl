struct NoVerticalDiffusion{NF} <: VerticalDiffusion{NF} end
NoVerticalDiffusion(SG::SpectralGrid) = NoVerticalDiffusion{SG.NF}()

function initialize!(   scheme::NoVerticalDiffusion,
                        model::PrimitiveEquation)
    return nothing
end 

# function barrier
function static_energy_diffusion!(  column::ColumnVariables,
                                    model::PrimitiveEquation)
    static_energy_diffusion!(column,model.static_energy_diffusion)
end

function static_energy_diffusion!(  column::ColumnVariables,
                                    scheme::NoVerticalDiffusion)
    return nothing
end

"""
Diffusion of dry static energy: A relaxation towards a reference
gradient of static energy wrt to geopotential, see Fortran SPEEDY documentation.
$(TYPEDFIELDS)"""
Base.@kwdef struct StaticEnergyDiffusion{NF<:AbstractFloat} <: VerticalDiffusion{NF}
    "time scale for strength"
    time_scale::Second = Hour(6)

    "[1] ∂SE/∂Φ, vertical gradient of static energy SE with geopotential Φ"
    static_energy_lapse_rate::NF = 0.1          
    
    # precomputations
    Fstar::Base.RefValue{NF} = Ref(zero(NF))    # excluding the surface pressure pₛ
end

StaticEnergyDiffusion(SG::SpectralGrid;kwargs...) = StaticEnergyDiffusion{SG.NF}(;kwargs...)

"""$(TYPEDSIGNATURES)
Initialize dry static energy diffusion."""
function initialize!(   scheme::StaticEnergyDiffusion,
                        model::PrimitiveEquation)

    (;nlev) = model.spectral_grid
    (;gravity) = model.planet
    C₀ = 1/nlev                     # average Δσ
    
    # Fortran SPEEDY documentation equation (70), excluding the surface pressure pₛ
    scheme.Fstar[] = C₀/gravity/scheme.time_scale.value
end

"""$(TYPEDSIGNATURES)
Apply dry static energy diffusion."""
function static_energy_diffusion!(  column::ColumnVariables,
                                    scheme::StaticEnergyDiffusion)
    
    (;nlev, dry_static_energy, flux_temp_upward, geopot) = column
    pₛ = column.pres[end]               # surface pressure
    Fstar = scheme.Fstar[]*pₛ
    Γˢᵉ = scheme.static_energy_lapse_rate 

    # relax static energy profile back to a reference gradient Γˢᵉ
    @inbounds for k in 1:nlev-1
        # Fortran SPEEDY doc eq (74)
        SEstar = dry_static_energy[k+1] + Γˢᵉ*(geopot[k] - geopot[k+1])
        SE = dry_static_energy[k]
        if SE < SEstar
            # SPEEDY code applies the flux to all layers below
            # effectively fluxing temperature from below the surface to
            # given height. It doesn't really make sense to me in
            # a meteorology sense to me since it can flux through stable layers,
            # but it's more stable and that's good enough for now.
            for kk = k+1:nlev
                flux_temp_upward[kk] += Fstar*(SEstar - SE)
            end
        end
    end
end

"""
Diffusion of dry static energy: A relaxation towards a reference
gradient of static energy wrt to geopotential, see Fortran SPEEDY documentation.
$(TYPEDFIELDS)"""
Base.@kwdef struct HumidityDiffusion{NF<:AbstractFloat} <: VerticalDiffusion{NF}
    "time scale for strength"
    time_scale::Second = Hour(24)

    "[1] ∂RH/∂σ, vertical gradient of relative humidity RH wrt sigma coordinate σ"
    humidity_gradient::NF = 0.5

    # precomputations
    Fstar::Base.RefValue{NF} = Ref(zero(NF))    # excluding the surface pressure pₛ
end

HumidityDiffusion(SG::SpectralGrid;kwargs...) = HumidityDiffusion{SG.NF}(;kwargs...)

"""$(TYPEDSIGNATURES)
Initialize dry static energy diffusion."""
function initialize!(   scheme::HumidityDiffusion,
                        model::PrimitiveEquation)

    (;nlev) = model.spectral_grid
    (;gravity) = model.planet
    C₀ = 1/nlev                     # average Δσ
    
    # Fortran SPEEDY documentation equation (70), excluding the surface pressure pₛ
    scheme.Fstar[] = C₀/gravity/scheme.time_scale.value

    # SPEEDY code version
    # scheme.Fstar[] = C₀/scheme.time_scale.value
end

# function barrier for all VerticalDiffusion, dispatch by type of humidity diffusion
function humidity_diffusion!(   column::ColumnVariables,
                                model::PrimitiveWet)
    humidity_diffusion!(column,model.humidity_diffusion,model)
end

# do nothing for primitive dry
function humidity_diffusion!(   column::ColumnVariables,
                                model::PrimitiveDry)
    return nothing
end

# do nothing for no vertical diffusion
function humidity_diffusion!(   column::ColumnVariables,
                                scheme::NoVerticalDiffusion,
                                model::PrimitiveEquation)
    return nothing
end

# function barrier to unpack model
function humidity_diffusion!(   column::ColumnVariables,
                                scheme::HumidityDiffusion,
                                model::PrimitiveEquation)
    humidity_diffusion!(column, scheme, model.geometry)
end

"""$(TYPEDSIGNATURES)
Apply humidity diffusion."""
function humidity_diffusion!(   column::ColumnVariables,
                                scheme::HumidityDiffusion,
                                geometry::Geometry)
    
    (;nlev, sat_humid, rel_humid, flux_humid_upward) = column
    pₛ = column.pres[end]               # surface pressure
    Fstar = scheme.Fstar[]*pₛ           
    # Fstar = scheme.Fstar[]              # SPEEDY code version
    Γ = scheme.humidity_gradient        # relative humidity wrt σ
    σ = geometry.σ_levels_full
    σ_half = geometry.σ_levels_half

    # relax humidity profile back to a reference gradient Γ
    # SPEEDY skips the uppermost levels, but maybe less relevant due to the
    # σ_half[k] scaling of the flux
    for k in 1:nlev-1  

        Γσ = Γ*(σ[k+1] - σ[k])
        Δrel_humid = rel_humid[k+1] - rel_humid[k]

        # Fortran SPEEDY doc eq (71)
        if Δrel_humid > Γσ
            # Fortran SPEEDY doc eq (72)
            # the σ_half[k] is not in the documentation, but it makes the diffusion
            # weaker higher up, so it's maybe not bad after all
            flux_humid_upward[k+1] += Fstar*σ_half[k]*sat_humid[k]*Δrel_humid
        
            # SPEEDY code version
            # flux = Fstar*σ_half[k]*sat_humid[k]*Δrel_humid
            # humid_tend[k] += flux*rsig[k]
            # humid_tend[k+1] -= flux*rsig[k+1]
        end
    end
end