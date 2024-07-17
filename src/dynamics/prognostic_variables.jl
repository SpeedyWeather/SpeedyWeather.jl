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
    (; trunc, nlat_half, nlev, nparticles) = SG
    (; NF, ArrayType) = SG
    (; SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector) = SG

    return PrognosticVariables{NF, ArrayType, nsteps,
        SpectralVariable2D, SpectralVariable3D, GridVariable2D, ParticleVector}(;
            trunc, nlat_half, nlayers=nlev, nparticles,
        )
end

"""$(TYPEDSIGNATURES)
Generator function."""
function PrognosticVariables(SG::SpectralGrid, model::ModelSetup)
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
    # dynamics
    progn_new.vor .= progn_old.vor
    progn_new.div .= progn_old.div
    progn_new.temp .= progn_old.temp
    progn_new.humid .= progn_old.humid
    progn_new.pres .= progn_old.pres

    # ocean
    progn_new.ocean.sea_surface_temperature .= progn_old.ocean.sea_surface_temperature
    progn_new.ocean.sea_ice_concentration .= progn_old.ocean.sea_ice_concentration
    
    # land
    progn_new.land.land_surface_temperature .= progn_old.land.land_surface_temperature
    progn_new.land.snow_depth .= progn_old.land.snow_depth
    progn_new.land.soil_moisture_layer1 .= progn_old.land.soil_moisture_layer1
    progn_new.land.soil_moisture_layer2 .= progn_old.land.soil_moisture_layer2

    progn_new.particles .= progn_old.particles
    progn_new.clock.time = progn_old.clock.time
    progn_new.scale[] = progn_old.scale[]

    return progn_new
end

"""
$(TYPEDSIGNATURES)
Sets new values for the keyword arguments (velocities, vorticity, divergence, etc..) into the
prognostic variable struct `progn` at timestep index `lf`. If `add==true` they are added to the 
current value instead. If a `SpectralTransform` S is provided, it is used when needed to set 
the variable, otherwise it is recomputed. In case `u` and `v` are set, `coslat_scaling_included`
specficies whether or not the 1/cos(lat) scaling is already included in the arrays or not (default:
`false`)

The input may be:
* A function or callable object `f(lond, latd, σ) -> value` or `f(lond, latd) -> value` (surface level variables)
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
    S::Union{Nothing, SpectralTransform} = nothing,
    coslat_scaling_included::Bool = false,
)
    isnothing(vor)   || set!(progn.vor[lf],     vor, geometry, S; add)
    isnothing(div)   || set!(progn.div[lf],     div, geometry, S; add)
    isnothing(temp)  || set!(progn.temp[lf],   temp, geometry, S; add)
    isnothing(humid) || set!(progn.humid[lf], humid, geometry, S; add)
    isnothing(pres)  || set!(progn.pres[lf],   pres, geometry, S; add)

    isnothing(sea_surface_temperature)  || set!(progn.ocean.sea_surface_temperature, sea_surface_temperature, geometry, S; add)
    isnothing(sea_ice_concentration)    || set!(progn.ocean.sea_ice_concentration, sea_ice_concentration, geometry, S; add)

    isnothing(land_surface_temperature) || set!(progn.land.land_surface_temperature, land_surface_temperature, geometry, S; add)
    isnothing(snow_depth)               || set!(progn.land.snow_depth, snow_depth, geometry, S; add)
    isnothing(soil_moisture_layer1)     || set!(progn.land.soil_moisture_layer1, soil_moisture_layer1, geometry, S; add)
    isnothing(soil_moisture_layer2)     || set!(progn.land.soil_moisture_layer2, soil_moisture_layer2, geometry, S; add)

    isnothing(u) | isnothing(v) || set_vordiv!(progn.vor[lf], progn.div[lf], u, v, geometry, S; add, coslat_scaling_included)
end

# set LTA <- LTA 
function set!(var::LowerTriangularArray{T}, L::LowerTriangularArray, varargs...; add::Bool) where T
    if add 
        if size(var) == size(L)
            var .+= T.(L) 
        else 
            L_var = spectral_truncation(L, size(var, 1, as=Matrix), size(var, 2, as=Matrix))
            var .+= L_var
        end 
    else 
        copyto!(var, L)
    end 
end 

# set LTA <- Grid 
function set!(var::LowerTriangularArray, grids::AbstractGridArray, geometry::Union{Geometry, Nothing}=nothing, S::Union{Nothing, SpectralTransform}=nothing; add)
    specs = isnothing(S) ? transform(grids) : transform(grids, S)
    set!(var, specs; add)
end

# set LTA <- func 
function set!(var::LowerTriangularArray, f::Function, geometry::Geometry{NF}, S::Union{SpectralTransform, Nothing}=nothing; add::Bool) where NF
    grid = zeros(geometry.Grid{NF}, geometry.nlat_half, geometry.nlev)
    set!(grid, f, geometry, S; add)

    if isnothing(S)
        spec = transform(grid)
        copyto!(var, spec)
    else 
        transform!(var, grid, S)
    end 
end

# set LTA <- number (change it to directly only set the [1,1] element, but not sure about the normalization)
function set!(var::LowerTriangularArray{T}, s::Number, geometry::Geometry{NF}, S::Union{SpectralTransform, Nothing}=nothing; add::Bool) where {T, NF}
    grid = ones(geometry.Grid{NF}, geometry.nlat_half, geometry.nlev) .* s
    s_spec = isnothing(S) ? transform(grid) : transform(grid, S)
    set!(var, s_spec, geometry, S; add)
end 

# set Grid <- Grid
function set!(var::AbstractGridArray, grids::AbstractGridArray, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add)
    if add 
        if grids_match(var, grids)
            var .+= grids
        else 
            var .+= interpolate(typeof(var), geometry.nlat_half, grids)
        end
    else 
        interpolate!(var, grids)
    end 
end 

# set Grid <- LTA
function set!(var::AbstractGridArray, specs::LowerTriangularArray, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add)
    grids = isnothing(S) ? transform(specs) : transform(specs, S)
    set!(var, grids, geometry, S; add)
end

# set Grid <- Func
function set!(var::AbstractGridArray, f::Function, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add)
    (; londs, latds, σ_levels_full) = geometry
    kernel(a, b) = add ? a+b : b
    for k in eachgrid(var)
        for ij in eachgridpoint(var)
            var[ij, k] = kernel(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
        end
    end
end

# set Grid (surface/single level) <- Func
function set!(var::AbstractGridArray{T,1}, f::Function, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add) where T
    (; londs, latds) = geometry
    kernel(a, b) = add ? a+b : b
    for ij in eachgridpoint(var)
        var[ij] = kernel(var[ij], f(londs[ij], latds[ij]))
    end
end

# set Grid <- Number 
function set!(var::AbstractGridArray, s::Number, geometry::Union{Geometry, Nothing}=nothing, S::Union{Nothing, SpectralTransform}=nothing; add::Bool)
    kernel(a, b) = add ? a+b : b
    sT = T(s)
    for k in eachgrid(var)
        for ij in eachgridpoint(var)
            var[ij, k] = kernel(var[ij, k], sT)
        end
    end
end 

# set Grid (surface/single level) <- Number 
function set!(var::AbstractGridArray{T,1}, s::Number, geometry::Union{Geometry, Nothing}=nothing, S::Union{Nothing, SpectralTransform}=nothing; add::Bool) where T
    kernel(a, b) = add ? a+b : b
    sT = T(s)
    for ij in eachgridpoint(var)
        var[ij] = kernel(var[ij], sT)
    end
end 

# set vor_div <- func 
function set_vordiv!(vor::LowerTriangularArray, div::LowerTriangularArray, u, v, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add::Bool, coslat_scaling_included::Bool=false) 
    
    u_L = similar(vor) 
    set!(u_L, u, geometry, S)
    v_L = similar(vor)
    set!(v_L, v, geometry, S)

    set_vordiv!(vor, div, u_L, v_L, geometry, S; add, coslat_scaling_included)
end

# set vor_div <- grid 
function set_vordiv!(vor::LowerTriangularArray, div::LowerTriangularArray, u::AbstractGridArray, v::AbstractGridArray, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add::Bool, coslat_scaling_included::Bool=false)
    
    if !coslat_scaling_included
        u_ = RingGrids.scale_coslat⁻¹(u)
        v_ = RingGrids.scale_coslat⁻¹(v)
    else
        u_ = u 
        v_ = v
    end 
    
    u_spec = isnothing(S) ? transform(u_) : transform(u_, S)
    v_spec = isnothing(S) ? transform(v_) : transform(v_, S)
    set_vordiv!(vor, div, u_spec, v_spec, geometry, S; add)
end 

# set vor_div <- LTA
function set_vordiv!(vor::LowerTriangularArray, div::LowerTriangularArray, u::LowerTriangularArray, v::LowerTriangularArray, geometry::Geometry, S::Union{Nothing, SpectralTransform}=nothing; add::Bool, coslat_scaling_included::Bool=false) 
  
    S = isnothing(S) ? SpectralTransform(geometry.spectral_grid) : S
    
    if !coslat_scaling_included
        u_ = transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)
        v_ = transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)
    else
        u_ = u 
        v_ = v
    end 

    if size(vor) != size(u_) != size(v_)
        u_new = zero(vor)
        copyto!(u_new, u_) 

        v_new = zero(vor)
        copyto!(v_new, v_)

        curl!(vor, u_new, v_new, S; add)
        divergence!(div, u_new, v_new, S; add)
    else 
        curl!(vor, u_, v_, S; add)
        divergence!(div, u_, v_, S; add)
    end
end 

set!(S::AbstractSimulation; kwargs...) = set!(S.prognostic_variables, S.model.geometry; S=S.model.spectral_transform, kwargs...)
