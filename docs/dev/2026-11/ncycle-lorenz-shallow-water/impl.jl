# experimental implicit formulations for NCycleLorenz + ShallowWater
using SpeedyWeather
const SW = SpeedyWeather
Base.@kwdef struct ImplicitSWX{NF} <: SW.AbstractImplicit
    centering::NF = 0.5
    μ::NF = 0.0          # 0: linear terms outside Lorenz average ("docs"), 1: Lawson (memory rotated)
end
SW.initialize!(::ImplicitSWX, ::SW.Variables, ::SW.AbstractModel) = nothing
SW.reinitialize!(::ImplicitSWX, args...) = nothing

function SW.implicit_correction!(vars::SW.Variables, implicit::ImplicitSWX, ts::SW.AbstractNCycleLorenz, model::SW.ShallowWater)
    dt = ts.Δt / vars.prognostic.scale[]
    w = SW.weight_coefficient(ts, vars.prognostic.clock)
    F_div = get_step(vars.tendencies.divergence, 1); G_div = get_step(vars.tendencies.divergence, 2)
    F_η = get_step(vars.tendencies.η, 1); G_η = get_step(vars.tendencies.η, 2)
    div = get_step(vars.prognostic.divergence, 1); η = get_step(vars.prognostic.η, 1)
    H = model.atmosphere.layer_thickness; g = model.planet.gravity
    α = implicit.centering; μ = implicit.μ
    ξ = α*dt
    l_indices = F_div.spectrum.l_indices
    @inbounds for lm in 1:size(F_div, 1)
        l = l_indices[lm]; ∇² = -l*(l-1)
        Sinv = inv(1 - ξ^2*H*g*∇²)
        solve(ad, aη) = (bd = Sinv*(ad - ξ*g*∇²*aη); (bd, aη - ξ*H*bd))
        xd = div[lm, 1]; xη = η[lm]
        Ld = -g*∇²*xη; Lη = -H*xd
        # new explicit memory
        Gd = w*(F_div[lm, 1] - Ld) + (1-w)*G_div[lm, 1]
        Gη = w*(F_η[lm] - Lη) + (1-w)*G_η[lm]
        # R x = x + dt S(Lx)
        sd, sη = solve(Ld, Lη); Rxd = xd + dt*sd; Rxη = xη + dt*sη
        # S G̃ and R G̃
        SGd, SGη = solve(Gd, Gη)
        tLd, tLη = -g*∇²*Gη, -H*Gd
        uLd, uLη = solve(tLd, tLη)
        RGd = Gd + dt*uLd; RGη = Gη + dt*uLη
        # x_new = R x + dt (μ R + (1-μ) S) G̃ ; G_after = (μ R + (1-μ) I) G̃
        xnd = Rxd + dt*(μ*RGd + (1-μ)*SGd); xnη = Rxη + dt*(μ*RGη + (1-μ)*SGη)
        Gad = μ*RGd + (1-μ)*Gd; Gaη = μ*RGη + (1-μ)*Gη
        # fold into F so update_prognostic gives G = G_after, x += dt*G_after
        F_div[lm, 1] = (Gad - (1-w)*G_div[lm, 1]) / w
        F_η[lm] = (Gaη - (1-w)*G_η[lm]) / w
        div[lm, 1] = xnd - dt*Gad
        η[lm] = xnη - dt*Gaη
    end
    return nothing
end

# existing formulation (linear terms inside Lorenz average) but per-substep centering θ_k
struct ImplicitSWθ{NF} <: SW.AbstractImplicit
    θ::Vector{NF}
end
SW.initialize!(::ImplicitSWθ, ::SW.Variables, ::SW.AbstractModel) = nothing
SW.reinitialize!(::ImplicitSWθ, args...) = nothing
function SW.implicit_correction!(vars::SW.Variables, implicit::ImplicitSWθ, ts::SW.AbstractNCycleLorenz, model::SW.ShallowWater)
    k = mod(vars.prognostic.clock.step_counter, length(implicit.θ)) + 1
    tmp = SW.ImplicitShallowWater{Float32}(implicit.θ[k])
    SW.implicit_correction!(vars, tmp, ts, model)
end

# experimental: implicit diffusion with Lorenz-weighted time step wΔt
const WDIFF = Ref(false)
function _wdiff!(tend, var, expl, w, Δt)
    l_indices = var.spectrum.l_indices
    @inbounds for lm in 1:size(tend, 1)
        l = l_indices[lm]; e = expl[l, 1]
        tend[lm, 1] = (tend[lm, 1] + e*var[lm, 1]) / (1 - w*Δt*e)
    end
end
function SW.horizontal_diffusion!(vars::SW.Variables, diffusion::SW.HyperDiffusion, model::SW.ShallowWaterModel)
    (WDIFF[] && model.time_stepping isa SW.NCycleLorenz) || return invoke(SW.horizontal_diffusion!, Tuple{SW.Variables, SW.AbstractHorizontalDiffusion, SW.ShallowWater}, vars, diffusion, model)
    ts = model.time_stepping
    w = SW.weight_coefficient(ts, vars.prognostic.clock)
    Δt = ts.Δt / model.planet.radius
    _wdiff!(get_step(vars.tendencies.vorticity, 1), get_step(vars.prognostic.vorticity, 1), diffusion.expl, w, Δt)
    _wdiff!(get_step(vars.tendencies.divergence, 1), get_step(vars.prognostic.divergence, 1), diffusion.expl_div, w, Δt)
    SW.zero_last_degree!(get_step(vars.tendencies.η, 1))
end
