function Base.show(io::IO, A::AbstractVariables)
    println(io, "$(typeof(A))")
    keys = propertynames(A)
    print_fields(io, A, keys)
end

# to be removed
struct PrognosticLayerTimesteps end
struct PrognosticSurfaceTimesteps end
struct PrognosticVariablesLayer end

export PrognosticVariablesOcean
@kwdef struct PrognosticVariablesOcean{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # OCEAN VARIABLES
    "Sea surface temperature [K]"
    sea_surface_temperature::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Sea ice concentration [1]"
    sea_ice_concentration::GridVariable2D = zeros(GridVariable2D, nlat_half)
end

export PrognosticVariablesLand
@kwdef struct PrognosticVariablesLand{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    GridVariable2D,         # <: AbstractGridArray
} <: AbstractPrognosticVariables
    # DIMENSION
    "Number of latitude rings on one hemisphere (Equator incl.), resolution parameter of grid"
    nlat_half::Int

    # LAND VARIABLES
    "Land surface temperature [K]"
    land_surface_temperature::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Snow depth [m]"
    snow_depth::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Soil moisture layer 1, volume fraction [1]"
    soil_moisture_layer1::GridVariable2D = zeros(GridVariable2D, nlat_half)

    "Soil moisture layer 2, volume fraction [1]"
    soil_moisture_layer2::GridVariable2D = zeros(GridVariable2D, nlat_half)
end

export PrognosticVariables
@kwdef struct PrognosticVariables{
    NF,                     # <: AbstractFloat
    ArrayType,              # Array, CuArray, ...
    NSTEPS,                 # number of timesteps
    SpectralVariable2D,     # <: LowerTriangularArray
    SpectralVariable3D,     # <: LowerTriangularArray
    GridVariable2D,         # <: AbstractGridArray
    ParticleVector,         # <: AbstractVector{Particle{NF}}
} <: AbstractPrognosticVariables

    # DIMENSIONS
    "max degree of spherical harmonics (0-based)"
    trunc::Int

    "Number of latitude rings on one hemisphere (Equator excl.), resolution parameter of grids"
    nlat_half::Int

    "number of vertical layers"
    nlayers::Int

    "Number of particles for particle advection"
    nparticles::Int

    # LAYERED VARIABLES
    "Vorticity of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    vor::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Divergence of horizontal wind field [1/s], but scaled by scale (=radius during simulation)"
    div::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Absolute temperature [K]"
    temp::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Specific humidity [kg/kg]"
    humid::NTuple{NSTEPS, SpectralVariable3D} =
        ntuple(i -> zeros(SpectralVariable3D, trunc+2, trunc+1, nlayers), NSTEPS)

    "Logarithm of surface pressure [log(Pa)] for PrimitiveEquation, interface displacement [m] for ShallowWaterModel"
    pres::NTuple{NSTEPS, SpectralVariable2D} =
        ntuple(i -> zeros(SpectralVariable2D, trunc+2, trunc+1), NSTEPS)

    "Random pattern following a random process [1]"
    random_pattern::SpectralVariable2D = zeros(SpectralVariable2D, trunc+2, trunc+1)

    "Ocean variables, sea surface temperature and sea ice concentration"
    ocean::PrognosticVariablesOcean{NF, ArrayType, GridVariable2D} =
        PrognosticVariablesOcean{NF, ArrayType, GridVariable2D}(; nlat_half)
    
    "Land variables, land surface temperature, snow and soil moisture"
    land::PrognosticVariablesLand{NF, ArrayType, GridVariable2D} =
        PrognosticVariablesLand{NF, ArrayType, GridVariable2D}(; nlat_half)

    "Particles for particle advection"
    particles::ParticleVector = zeros(ParticleVector, nparticles)

    "Scaling for vor, div. scale=1 outside simulation, =radius during simulation"
    scale::Base.RefValue{NF} = Ref(one(NF))

    "Clock that keeps track of time, number of timesteps to integrate for."
    clock::Clock = Clock()
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid; nsteps=DEFAULT_NSTEPS)
    (; trunc, nlat_half, nlayers, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector) = SG

    return PrognosticVariables{NF, ArrayType, nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector}(;
            trunc, nlat_half, nlayers=nlayers, nparticles,
        )
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid, model::AbstractModel)
    PrognosticVariables(SG, nsteps = model.time_stepping.nsteps)
end

function Base.show(
    io::IO,
    progn::PrognosticVariables{NF, ArrayType, NSTEPS},
) where {NF, ArrayType, NSTEPS}
    Grid = typeof(progn.ocean.sea_surface_temperature)
    println(io, "PrognosticVariables{$NF, $ArrayType}")
    
    # variables
    (; trunc, nlat_half, nlayers, nparticles) = progn
    nlat = RingGrids.get_nlat(Grid, nlat_half)
    println(io, "├ vor:   T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ div:   T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ temp:  T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ humid: T$trunc, $nlayers-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ pres:  T$trunc, 1-layer, $NSTEPS-steps LowerTriangularArray{$NF}")
    println(io, "├ random_pattern: T$trunc, 1-layer LowerTriangularArray{$NF}")
    println(io, "├┐ocean: PrognosticVariablesOcean{$NF}")
    println(io, "│├ sea_surface_temperature:  $nlat-ring $Grid")
    println(io, "│└ sea_ice_concentration:    $nlat-ring $Grid")
    println(io, "├┐land:  PrognosticVariablesLand{$NF}")
    println(io, "│├ land_surface_temperature: $nlat-ring $Grid")
    println(io, "│├ snow_depth: $nlat-ring $Grid")
    println(io, "│├ soil_moisture_layer1:     $nlat-ring $Grid")
    println(io, "│└ soil_moisture_layer2:     $nlat-ring $Grid")
    println(io, "├ particles: $nparticles-element $(typeof(progn.particles))")
    println(io, "├ scale: $(progn.scale[])")
    print(io,   "└ clock: $(progn.clock.time)")
end

# has(::PrognosticVariables{NF, Grid, M}, var_name::Symbol) where {NF, Grid, M} = has(M, var_name)

"""$(TYPEDSIGNATURES)
Copies entries of `progn_old` into `progn_new`."""
function Base.copy!(progn_new::PrognosticVariables, progn_old::PrognosticVariables)

    for i in eachindex(progn_new.vor)   # each leapfrog time step
        progn_new.vor[i] .= progn_old.vor[i]
        progn_new.div[i] .= progn_old.div[i]
        progn_new.temp[i] .= progn_old.temp[i]
        progn_new.humid[i] .= progn_old.humid[i]
        progn_new.pres[i] .= progn_old.pres[i]
    end

    # ocean
    progn_new.ocean.sea_surface_temperature .= progn_old.ocean.sea_surface_temperature
    progn_new.ocean.sea_ice_concentration .= progn_old.ocean.sea_ice_concentration
    
    # land
    progn_new.land.land_surface_temperature .= progn_old.land.land_surface_temperature
    progn_new.land.snow_depth .= progn_old.land.snow_depth
    progn_new.land.soil_moisture_layer1 .= progn_old.land.soil_moisture_layer1
    progn_new.land.soil_moisture_layer2 .= progn_old.land.soil_moisture_layer2

    # copy largest subset of particles
    if length(progn_new.particles) != length(progn_old.particles)
        nnew = length(progn_new.particles)
        nold = length(progn_old.particles)
        nsub = min(nnew, nold)
        @warn "Number of particles changed (origin: $nold, destination: $nnew), copying over only the largest subset ($nsub particles)"
        progn_new.particles[1:nsub] .= progn_old.particles[1:nsub]
    else
        progn_new.particles .= progn_old.particles
    end

    progn_new.clock.time = progn_old.clock.time
    progn_new.scale[] = progn_old.scale[]

    return progn_new
end

export set!

"""
$(TYPEDSIGNATURES)
Sets new values for the keyword arguments (velocities, vorticity, divergence, etc..) into the
prognostic variable struct `progn` at timestep index `lf`. If `add==true` they are added to the 
current value instead. If a `SpectralTransform` S is provided, it is used when needed to set 
the variable, otherwise it is recomputed. In case `u` and `v` are provied, actually the divergence
and vorticity are set and `coslat_scaling_included` specficies whether or not the 1/cos(lat) 
scaling is already included in the arrays or not (default: `false`)

The input may be:
* A function or callable object `f(lond, latd, σ) -> value` (multilevel variables) 
* A function or callable object `f(lond, latd) -> value` (surface level variables)
* An instance of `AbstractGridArray` 
* An instance of `LowerTriangularArray` 
* A scalar `<: Number` (interpreted as a constant field in grid space)
"""
function set!(
    progn::PrognosticVariables,
    geometry::Geometry;
    u = nothing,
    v = nothing,
    vor = nothing,
    div = nothing,
    temp = nothing,
    humid = nothing,
    pres = nothing,
    sea_surface_temperature = nothing, 
    sea_ice_concentration = nothing, 
    land_surface_temperature = nothing, 
    snow_depth = nothing, 
    soil_moisture_layer1 = nothing, 
    soil_moisture_layer2 = nothing,
    lf::Integer = 1,
    add::Bool = false,
    spectral_transform::Union{Nothing, SpectralTransform} = nothing,
    coslat_scaling_included::Bool = false,
)
    # ATMOSPHERE
    isnothing(vor)   || set!(progn.vor[lf],     vor, geometry, spectral_transform; add)
    isnothing(div)   || set!(progn.div[lf],     div, geometry, spectral_transform; add)
    isnothing(temp)  || set!(progn.temp[lf],   temp, geometry, spectral_transform; add)
    isnothing(humid) || set!(progn.humid[lf], humid, geometry, spectral_transform; add)
    isnothing(pres)  || set!(progn.pres[lf],   pres, geometry, spectral_transform; add)
    
    # or provide u, v instead of vor, div
    isnothing(u) | isnothing(v) || set_vordiv!(progn.vor[lf], progn.div[lf], u, v, geometry, spectral_transform; add, coslat_scaling_included)
    
    # OCEAN
    isnothing(sea_surface_temperature)  || set!(progn.ocean.sea_surface_temperature, sea_surface_temperature, geometry, spectral_transform; add)
    isnothing(sea_ice_concentration)    || set!(progn.ocean.sea_ice_concentration, sea_ice_concentration, geometry, spectral_transform; add)

    # LAND
    isnothing(land_surface_temperature) || set!(progn.land.land_surface_temperature, land_surface_temperature, geometry, spectral_transform; add)
    isnothing(snow_depth)               || set!(progn.land.snow_depth, snow_depth, geometry, spectral_transform; add)
    isnothing(soil_moisture_layer1)     || set!(progn.land.soil_moisture_layer1, soil_moisture_layer1, geometry, spectral_transform; add)
    isnothing(soil_moisture_layer2)     || set!(progn.land.soil_moisture_layer2, soil_moisture_layer2, geometry, spectral_transform; add)
    return nothing
end

# set LTA <- LTA 
function set!(
    var::LowerTriangularArray,
    L::LowerTriangularArray,
    varargs...;
    add::Bool=false,
)
    if add 
        if size(var) == size(L)
            var .+= L
        else 
            L_var = spectral_truncation(L, size(var, 1, as=Matrix), size(var, 2, as=Matrix))
            var .+= L_var
        end 
    else 
        size(var) != size(L) || fill!(var, 0) # copyto! copies over the largest subset, when size(var) > size(L), the copyto! isn't enough by itself
        copyto!(var, L)
    end 
    return var
end 

# set LTA <- Grid 
function set!(
    var::LowerTriangularArray,
    grids::AbstractGridArray,
    geometry::Union{Geometry, Nothing}=nothing,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
)
    if isnothing(S)
        specs = transform(grids)
    else
        # convert to number format in S, needed for FFTW
        grids = convert.(eltype(S), grids)
        specs = transform(grids, S)
    end
    set!(var, specs; add)
end

# set LTA <- func 
function set!(
    var::LowerTriangularArray,
    f::Function,
    geometry::Geometry{NF, Grid},
    S::Union{SpectralTransform, Nothing}=nothing;
    add::Bool=false,
) where {NF, Grid}
    grid = ndims(var) == 1 ? zeros(Grid{NF}, geometry.nlat_half) : zeros(Grid{NF}, geometry.nlat_half, geometry.nlayers)
    set!(grid, f, geometry, S; add=false)
    set!(var, grid, geometry, S; add)
end

# set LTA <- number
function set!(
    var::LowerTriangularArray{T},
    s::Number,
    geometry::Geometry{NF},
    S::Union{SpectralTransform, Nothing}=nothing;
    add::Bool=false,
) where {T, NF}
    
    # appropiate normalization, assume standard 2√π normalisation if no transform is given 
    norm_sphere = isnothing(S) ? 2sqrt(π) : S.norm_sphere

    # all elements are zero except for the 0,0 one
    var_new = zero(var)

    for k in eachmatrix(var_new)
        var_new[1, k] = norm_sphere * s
    end 

    set!(var, var_new, geometry, S; add)
end 

# set Grid <- Grid
function set!(
    var::AbstractGridArray,
    grids::AbstractGridArray,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
)
    if add 
        if grids_match(var, grids)
            var .+= grids
        else 
            var .+= interpolate(typeof(var), geometry.nlat_half, grids)
        end
    else 
        interpolate!(var, grids)
    end 
    return var 
end 

# set Grid <- LTA
function set!(
    var::AbstractGridArray,
    specs::LowerTriangularArray,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
)
    grids = isnothing(S) ? transform(specs) : transform(specs, S)
    set!(var, grids, geometry, S; add)
end

# set Grid <- Func
function set!(
    var::AbstractGridArray,
    f::Function,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
)
    (; londs, latds, σ_levels_full) = geometry
    kernel(a, b) = add ? a+b : b
    for k in eachgrid(var)
        for ij in eachgridpoint(var)
            var[ij, k] = kernel(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
        end
    end
    return var
end

# set Grid (surface/single level) <- Func
function set!(
    var::AbstractGridArray{T,1},
    f::Function,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
) where T
    (; londs, latds) = geometry
    kernel(a, b) = add ? a+b : b
    for ij in eachgridpoint(var)
        var[ij] = kernel(var[ij], f(londs[ij], latds[ij]))
    end
    return var
end

# set Grid <- Number 
function set!(
    var::AbstractGridArray{T}, 
    s::Number, 
    geometry::Union{Geometry, Nothing}=nothing, 
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
) where T
    kernel(a, b) = add ? a+b : b
    sT = T(s)
    var .= kernel.(var, sT)
end 

# set vor_div <- func 
function set_vordiv!(
    vor::LowerTriangularArray,
    div::LowerTriangularArray,
    u_func,
    v_func,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
    coslat_scaling_included::Bool=false,
)
    u_L = similar(vor) 
    set!(u_L, u_func, geometry, S)
    v_L = similar(vor)
    set!(v_L, v_func, geometry, S)

    set_vordiv!(vor, div, u_L, v_L, geometry, S; add, coslat_scaling_included)
end

# set vor_div <- grid 
function set_vordiv!(
    vor::LowerTriangularArray,
    div::LowerTriangularArray,
    u::AbstractGridArray,
    v::AbstractGridArray,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
    coslat_scaling_included::Bool=false,
)
    u_ = coslat_scaling_included ? u : RingGrids.scale_coslat⁻¹(u)
    v_ = coslat_scaling_included ? v : RingGrids.scale_coslat⁻¹(v)

    u_spec = isnothing(S) ? transform(u_) : transform(u_, S)
    v_spec = isnothing(S) ? transform(v_) : transform(v_, S)

    set_vordiv!(vor, div, u_spec, v_spec, geometry, S; add, coslat_scaling_included=true)
end 

# set vor_div <- LTA
function set_vordiv!(
    vor::LowerTriangularArray,
    div::LowerTriangularArray,
    u::LowerTriangularArray,
    v::LowerTriangularArray,
    geometry::Geometry,
    S::Union{Nothing, SpectralTransform}=nothing;
    add::Bool=false,
    coslat_scaling_included::Bool=false,
) 
    S = isnothing(S) ? SpectralTransform(geometry.spectral_grid) : S
     
    u_ = coslat_scaling_included ? u : transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)
    v_ = coslat_scaling_included ? v : transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)

    if size(vor) != size(u_) != size(v_)
        u_new = zero(vor)
        copyto!(u_new, u_) 

        v_new = zero(vor)
        copyto!(v_new, v_)

        curl!(vor, u_new, v_new, S; add, radius=geometry.radius)
        divergence!(div, u_new, v_new, S; add, radius=geometry.radius)
    else 
        curl!(vor, u_, v_, S; add, radius=geometry.radius)
        divergence!(div, u_, v_, S; add, radius=geometry.radius)
    end
end 

"""
$(TYPEDSIGNATURES)

Sets properties of the simuluation `S`. Convenience wrapper to call the other concrete 
`set!` methods. All `kwargs` are forwarded to these methods, which are documented 
seperately. See their documentation for possible `kwargs`. 
"""
function set!(S::AbstractSimulation; kwargs...)
    set!(S.prognostic_variables, S.model.geometry; spectral_transform=S.model.spectral_transform, kwargs...)
end

function set!(progn::PrognosticVariables, model::AbstractModel; kwargs...)
    progn.scale[] != 1 && @warn "Prognostic variables are scaled with $(progn.scale[]), but `set!` assumes unscaled variables."
    set!(progn, model.geometry; spectral_transform=model.spectral_transform, kwargs...)
end