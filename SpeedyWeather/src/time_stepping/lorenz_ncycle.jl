export NCycleLorenz, NCycleLorenzA, NCycleLorenzB, NCycleLorenzAB, NCycleLorenzABBA

"""Abstract type for Lorenz N-cycle variants following Hotta et al. (2016)."""
abstract type NCycleLorenzVariant end

"""Version A: weights w_k = N/(N-k) for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzA <: NCycleLorenzVariant end

"""Version B: weights w_k = N/k for k=1,...,N-1; w_0 = 1"""
struct NCycleLorenzB <: NCycleLorenzVariant end

"""Version AB: alternates A and B every N steps"""
struct NCycleLorenzAB <: NCycleLorenzVariant end

"""Version ABBA: uses A-B-B-A sequence (only for N=4, provides 4th order accuracy)"""
struct NCycleLorenzABBA <: NCycleLorenzVariant end

"""
    NCycleLorenz{NF, V} <: AbstractTimeStepper

A semi-implicit Lorenz N-cycle time integration scheme following Hotta et al. (2016).

# Algorithm (per substep)
1. G = w*F_E(x) + (1-w)*G     (weighted tendency accumulation)
2. dx = (I - α*Δt*L_I)^(-1) * (G + L_I*x)  (implicit solve)
3. x = x + Δt*dx              (state update)
$(TYPEDFIELDS)
"""
mutable struct NCycleLorenz{NF, V, IntType, S, MS, B} <: AbstractTimeStepper
    "[OPTION] Number of cycles N (3 or 4 recommended, 4 is more stable)"
    cycles::IntType
    
    "[OPTION] Variant: NCycleLorenzA() (default), B, AB, or ABBA"
    variant::V
    
    "[OPTION] Time step for T31, scale linearly with resolution"
    Δt_at_T31::S

    "[OPTION] Adjust `Δt_at_T31` with the `output_dt` to reach `output_dt` exactly"
    adjust_with_output::B

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt_sec::NF

    "[DERIVED] Time step Δt [s/m] at specified resolution, scaled by 1/radius"
    Δt::NF
end

function Adapt.adapt_structure(to, L::NCycleLorenz)
    return (; Δt = L.Δt, Δt_sec = L.Δt_sec, Δt_millisec = L.Δt_millisec)
end

# first for prognostic variables, 2nd for "prognostic" tendency from previous time step
get_prognostic_steps(::NCycleLorenz) = 2

"""$(TYPEDSIGNATURES)
Generator function for NCycleLorenz struct using `spectral_grid` for resolution."""
function NCycleLorenz(
    spectral_grid::SpectralGrid;
    cycles = 4,
    variant = NCycleLorenzA(),
    Δt_at_T31 = Minute(30),
    adjust_with_output = true,
    radius = DEFAULT_RADIUS,
)
    (; NF, trunc) = spectral_grid
    
    # compute time step
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, DEFAULT_RADIUS, adjust_with_output)
    Δt_sec::NF = Δt_millisec.value/1000
    Δt::NF = Δt_sec/radius

    return NCycleLorenz(cycles, variant, Second(Δt_at_T31), adjust_with_output, Δt_millisec, Δt_sec, Δt)
end

"""$(TYPEDSIGNATURES)
Get current substep within N-cycle (0 to N-1) from the clock."""
current_substep(L::NCycleLorenz, clock) = mod(clock.timestep_counter, L.cycles)

"""$(TYPEDSIGNATURES)
Initialize NCycleLorenz time stepper."""
function initialize!(L::NCycleLorenz, model::AbstractModel)
    if L.variant isa NCycleLorenzABBA && L.cycles != 4          # Validate compatibility
        @warn "N-Cycle Lorenz with ABBA variant is for N=4 (4th order accurate), but N=$(L.cycles). Consider cycles=4 or variant A/B/AB."
    end

    calculate_timestep!(L, model)
    return nothing
end

# Weight coefficient functions (Hotta et al. 2016, Eqs 2-3, 7-8)
"""$(TYPEDSIGNATURES)
Compute weight coefficient w for current substep."""
function weight_coefficient(L::NCycleLorenz{NF}, clock) where NF
    k = current_substep(L, clock)
    return weight_coefficient(NF, L.variant, k, L.cycles)
end

# Type-stable versions with explicit NF
@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzA, k::Int, N::Int) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, N - k)
end

@inline function weight_coefficient(::Type{NF}, ::NCycleLorenzB, k::Int, N::Int) where NF
    return k == 0 ? one(NF) : convert(NF, N) / convert(NF, k)
end

# Fallback for when NF is not specified (uses DEFAULT_NF)
@inline weight_coefficient(V::NCycleLorenzVariant, k::Int, N::Int) = 
    weight_coefficient(DEFAULT_NF, V, k, N)

"""$(TYPEDSIGNATURES)
Accumulate weighted tendency: G = w*F_E + (1-w)*G (Hotta et al. 2016, Eq. 19)"""
function accumulate_tendency!(
    G::LowerTriangularArray,
    F_explicit::LowerTriangularArray,
    w::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(G) == size(F_explicit) || throw(BoundsError())
    
    w_NF = convert(NF, w)
    
    launch!(
        architecture(G), 
        SpectralWorkOrder, 
        size(G), 
        accumulate_tendency_kernel!, 
        G, F_explicit, w_NF
    )
    
    return nothing
end

@kernel inbounds=true function accumulate_tendency_kernel!(
    G, F_explicit, @Const(w)
)
    lmk = @index(Global, Linear)
    G[lmk] = w * F_explicit[lmk] + (1 - w) * G[lmk]
end

"""$(TYPEDSIGNATURES)
Update state: x = x + Δt*G (Hotta et al. 2016, Eq. 21)"""
function update_state!(
    x::LowerTriangularArray,
    G::LowerTriangularArray,
    dt::Real,
    L::NCycleLorenz{NF},
) where NF
    @boundscheck size(x) == size(G) || throw(BoundsError())
    
    dt_NF = convert(NF, dt)
    
    launch!(
        architecture(x), 
        SpectralWorkOrder, 
        size(x), 
        update_state_kernel!, 
        x, G, dt_NF
    )
    
    return nothing
end

@kernel inbounds=true function update_state_kernel!(
    x, G, @Const(dt)
)
    lmk = @index(Global, Linear)
    x[lmk] = x[lmk] + dt * G[lmk]
end

"""$(TYPEDSIGNATURES)
Lorenz N-cycle time stepping for all prognostic variables following Hotta et al. (2016).

This function performs ONE substep of the N-cycle. The clock tracks which substep we're on.

Algorithm per substep:
1. Accumulate weighted tendency: G = w*F_E + (1-w)*G (Eq. 19)
2. Apply implicit correction: G → (I - α*Δt*L_I)^(-1) * (G + L_I*x) (Eq. 20)
3. Update state: x = x + Δt*G (Eq. 21)
"""
function lorenz_ncycle_timestep!(
    vars::Variables,
    dt::Real,
    model::AbstractModel,
)
    L = model.time_stepping::NCycleLorenz
    clock = progn.clock
    
    # Get current weight from the clock
    w_current = weight_coefficient(L, clock)
    
    # Step 1: Accumulate weighted explicit tendencies (Eq. 19)
    for (varname, tendname) in zip(prognostic_variables(model), tendency_names(model))
        var = getfield(progn, varname)
        G = get_step(var, 2)
        
        F_explicit = getfield(tend, tendname)
        spectral_truncation!(F_explicit)
        
        accumulate_tendency!(G, F_explicit, w_current, L)
    end
    
    # Accumulate for tracers
    for (name, tracer) in model.tracers
        if tracer.active
            var = progn.tracers[name]
            G = get_step(var, 2)
            
            F_explicit = tend.tracers_tend[name]
            spectral_truncation!(F_explicit)
            
            accumulate_tendency!(G, F_explicit, w_current, L)
        end
    end
    
    # Step 2: Apply implicit correction (Eq. 20)
    lorenz_implicit_correction!(diagn, progn, model.implicit, model)

    # Step 3: Update states with corrected tendencies (Eq. 21)
    for varname in prognostic_variables(model)
        var = getfield(progn, varname)
        x = get_step(var, 1)
        G = get_step(var, 2)
        
        update_state!(x, G, dt, L)
    end
    
    # Update tracers
    for (name, tracer) in model.tracers
        if tracer.active
            var = progn.tracers[name]
            x = get_step(var, 1)
            G = get_step(var, 2)
            
            update_state!(x, G, dt, L)
        end
    end
    
    
    # Evolve random pattern
    random_process!(progn, model.random_process)
    
    return nothing
end