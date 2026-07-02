export Leapfrog

abstract type AbstractLeapfrog <: AbstractTimeStepper end

"""Leapfrog time stepping defined by the following fields
$(TYPEDFIELDS)"""
mutable struct Leapfrog{NF, S, B, MS} <: AbstractLeapfrog
    "[OPTION] Time step for T31, scale linearly to spectral resolution `trunc`"
    Δt_at_T31::S

    "[OPTION] Adjust `Δt_at_T31` with the output `interval` to output exactly after integer time steps"
    adjust_with_output::B

    "[OPTION] Robert (1966) time filter coefficient to suppress the computational mode"
    robert_filter::NF

    "[OPTION] Williams time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::NF

    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt::NF
end

Adapt.adapt_structure(to, L::Leapfrog) = Adapt.adapt_structure(to, LeapfrogCore(L.Δt_millisec, L.Δt))

# HOW MANY STEPS DO VARIABLES NEED?
# leapfrogging always needs 2 steps in spectral
prognostic_spectral_steps(::AbstractLeapfrog) = 2
# but in 2D only 1 step in grid space
prognostic_grid_steps(::AbstractLeapfrog, ::Union{<:Barotropic, <:ShallowWater}) = 1
# but the parameterizations are evaluated at the previous step so 2
prognostic_grid_steps(::AbstractLeapfrog, ::PrimitiveEquation) = 2
# always only one step for tendencies
tendency_steps(::AbstractLeapfrog) = 1

# WHICH STEP TO READ WHEN
# in Leapfrog use the current (=2nd) step to do the transforms
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::SpeedyTransforms.AbstractSpectralTransform) = 2

# but in Barotropc/ShallowWater models the 2nd one doesn't exist on the grid and the 1st is considered to be the current step
@inline which_prognostic_step(var::AbstractField, ::AbstractLeapfrog, ::SpeedyTransforms.AbstractSpectralTransform, ::TwoDModels) = 1
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractForcing) = 2
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractDrag) = 2

# in Leapfrog use the current (=2nd) in the dynamical core
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractDynamicalCoreComponent) = 2

# the linear terms in the dynamical core have to be evaluated at the previous time step
# they are then moved forward within the implicit corrections
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractGeopotential) = 1
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::LinearDynamicalCore) = 1

# but in Barotropc/ShallowWater models the 2nd one doesn't exist on the grid and the 1st is considered to be the current step
@inline which_prognostic_step(var::AbstractField, ::AbstractLeapfrog, ::AbstractDynamicalCoreComponent, ::TwoDModels) = 1

@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractHorizontalDiffusion) = 1
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::DiffusiveVerticalAdvection) = 1
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::DispersiveVerticalAdvection) = 2

# Parameterizations should always be evaluated on the previous time step for Euler forward
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractParameterization) = 1

# by default use the first step here but you may extend for subtypes of AbstractOcean elsewhere
# e.g. to leapfrog also a SlabOcean model's SST
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractOcean) = 1
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractSeaIce) = 1

# particle advection using u, v at current not previous time step
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractParticleAdvection) = 2

# dispatch over timestepper to decide between implicit or explicit diffusion
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::Nothing, ::AbstractLeapfrog) = true
@inline implicit_diffusion(::AbstractHorizontalDiffusion, ::AbstractImplicit, ::AbstractLeapfrog) = true

# particle advection reads u, v, use 1st and only step in 2D models, but 2nd and current step in 3D models
@inline which_prognostic_step(var::AbstractField, ::AbstractLeapfrog, ::AbstractParticleAdvection) = 2
@inline which_prognostic_step(var::AbstractField, ::AbstractLeapfrog, ::AbstractParticleAdvection, ::TwoDModels) = 1

# copy prognostic variables' 1st step to 2nd, that way which_prognostic_step can always be 2
function initialize!(vars::Variables, ::AbstractLeapfrog, ::AbstractModel)
    # time step variables are dynamically defined by existence in tendencies
    # but statically compiled into the tendency_names function
    (; prognostic) = vars
    for varname in tendency_names(vars)
        var_old, var_new = get_steps(getfield(prognostic, varname))
        var_new .= var_old
    end
    for varname in tracer_tendency_names(vars)
        var_old, var_new = get_steps(getfield(prognostic.tracers, varname))
        var_new .= var_old
    end
    for varname in ocean_tendency_names(vars)
        var_old, var_new = get_steps(getfield(prognostic.ocean, varname))
        var_new .= var_old
    end
    for varname in land_tendency_names(vars)
        var_old, var_new = get_steps(getfield(prognostic.land, varname))
        var_new .= var_old
    end
    return nothing
end

# copy grid prognostics from 2nd to 1st step to retain these fields for the previous time step
function move_prognostic_grid_variables_back!(
        vars::Variables,
        time_stepping::AbstractLeapfrog,
        model::AbstractModel,
    )
    # The atmospheric grid prognostics (vorticity, divergence, temperature, [humidity], pressure/η)
    # live in the `:grid` fuse parent and the velocities (u, v) in the `:uv_grid` fuse parent.
    # Together these cover exactly `tendency_and_uv_names`, so a single contiguous copy from the
    # current step (2nd) to the previous step (1st) per fuse parent replaces the per-variable loop.
    (; grid) = vars
    for fused in (vars.fused.grid, vars.fused.uv_grid)
        var_old, var_new = get_steps(parent(fused))
        var_old .= var_new
    end
    # tracers, ocean and land grid variables are not part of the atmospheric fuse — copy per-variable
    for varname in tracer_tendency_names(vars)
        var_old, var_new = get_steps(getfield(grid.tracers, varname))
        var_old .= var_new
    end
    # does not have to be done for ocean/land as they don't have spectral variables
    # the grid variables are already in prognostic and the time stepper takes care
    # of moving them back during the update
    return nothing
end

# for leapfrog do first semi-implicit corrections then horizontal diffusion
function diffusion_and_implicit!(vars, ::AbstractLeapfrog, ::AbstractImplicit, model)
    implicit_correction!(vars, model)
    horizontal_diffusion!(vars, model)
    return nothing
end

"""($TYPEDSIGNATURES)
Immutable core struct used to adapt only the time step fields for use in kernels"""
struct LeapfrogCore{NF, MS} <: AbstractLeapfrog
    "[DERIVED] Time step Δt in milliseconds at specified resolution"
    Δt_millisec::MS

    "[DERIVED] Time step Δt [s] at specified resolution"
    Δt::NF
end

Adapt.@adapt_structure LeapfrogCore

"""$(TYPEDSIGNATURES)
Generator function for a Leapfrog struct using `spectral_grid`
for the resolution information."""
function Leapfrog(
        spectral_grid::SpectralGrid;
        Δt_at_T31 = Minute(40),
        adjust_with_output = true,
        robert_filter = 0.1,
        williams_filter = 0.53,
    )
    (; NF, trunc) = spectral_grid

    # compute time step
    Δt_millisec::Millisecond = get_Δt_millisec(Second(Δt_at_T31), trunc, DEFAULT_RADIUS, adjust_with_output)
    Δt::NF = Δt_millisec.value / 1000

    return Leapfrog(
        Second(Δt_at_T31), adjust_with_output, NF(robert_filter), NF(williams_filter), Δt_millisec, Δt,
    )
end

"""$(TYPEDSIGNATURES)
Initialize leapfrogging `L` by recalculating the time step given the output time step
`interval` from `model.output`. Recalculating will slightly adjust the time step to
be a divisor such that an integer number of time steps matches exactly with the output
time step."""
function initialize!(L::Leapfrog, model::AbstractModel)
    calculate_Δt!(L, model)         # common among several time steppers
    return nothing
end

"""$(TYPEDSIGNATURES) Leapfrog is spun up with 1 Euler forward step that doesn't count for clock + output"""
spin_up_steps(::AbstractLeapfrog) = 1

function time_step!(clock::Clock, time_stepping::Leapfrog)
    Δt = time_stepping.Δt_millisec  # ::Millisecond, integer based hence ÷ not / below
    i = clock.step_counter          # 0-based as the clock is only stepped below
    @trace if i == 0                # first Euler step at Δt/2
        # i counts every time step, for the clock the first Euler step does not count
        # hence after this the time_stepping will be 1 ahead of clock step counter
        time_step!(clock, Δt ÷ 2, increase_counter = false)
    elseif i == 1                   # second step: Leapfrog at Δt
        # subtract the Δt/2 again as otherwise the time can be 1ms off due to rounding
        clock.time -= Δt ÷ 2
        time_step!(clock, Δt)
    else                            # later steps: Leapfrog at 2Δt but increase clock by Δt
        time_step!(clock, Δt)
    end
    return nothing
end

default_time_step(L::Leapfrog) = 2 * L.Δt
function time_step(L::Leapfrog, clock::Clock)
    (; Δt) = L
    return ifelse(
        clock.step_counter == 0, Δt / 2,        # first step Euler with Δt/2
        ifelse(
            clock.step_counter == 1, Δt,        # 2nd step leapfrog with Δt
            default_time_step(L)                # later steps leapfrog with 2Δt
        )
    )
end

# on 1st Euler step + 2nd Leapfrog step disable filters; with RAW filters afterwards
prognostic_step(::Leapfrog, clock::Clock) = ifelse(clock.step_counter <= 1, 1, 2)

function update_prognostic!(
        var::AbstractArray,
        tendency::AbstractArray,
        clock::Clock,
        time_stepping::Leapfrog,
        implicit::Union{Nothing, AbstractImplicit},
        ::AbstractModel,
        scale::Real = 1,
    )
    Δt = time_step(time_stepping, clock)
    Δt /= oftype(Δt, scale)                         # scale time step on the fly *1/radius for atmospheric variables
    lf = prognostic_step(time_stepping, clock)      # leapfrog prognostic step index
    var_old, var_new = get_steps(var)
    var_lf = get_step(var, lf)                      # view on either t or t+dt to dis/enable Williams filter
    var_tend = get_tendency_step(tendency, time_stepping, time_stepping)

    @boundscheck lf == 1 || lf == 2 || throw(BoundsError())
    @boundscheck size(var_old) == size(var_new) == size(var_tend) || throw(BoundsError())

    (; robert_filter, williams_filter) = time_stepping          # coefficients for filters

    # LEAP FROG time step with or without Robert+Williams filter
    # Robert time filter to compress computational mode, Williams filter for 3rd order accuracy
    # see Williams (2009), Eq. 7-9
    # for lf == 1 (time steps 1 or 2) no filters applied (w1=w2=0)
    # for lf == 2 (later steps) Robert+Williams filter is applied
    w1 = (lf - 1) * robert_filter * williams_filter / 2         # = ν*α/2 in Williams (2009, Eq. 8)
    w2 = (lf - 1) * robert_filter * (1 - williams_filter) / 2   # = ν(1-α)/2 in Williams (2009, Eq. 9)

    launch!(
        architecture(var_tend), SpectralWorkOrder, size(var_tend), leapfrog_kernel!,
        var_old, var_new, var_lf, var_tend, Δt, w1, w2
    )
    return nothing
end

@kernel inbounds = true function leapfrog_kernel!(var_old, var_new, var_lf, tendency, Δt, w1, w2)
    lmk = @index(Global, Linear)    # every harmonic lm, every vertical layer k
    old = var_old[lmk]
    new = old + Δt * tendency[lmk]
    update = old - 2var_lf[lmk] + new
    var_old[lmk] = var_lf[lmk] + w1 * update
    var_new[lmk] = new - w2 * update
end
