@testset "Callbacks interface" begin
    
    Base.@kwdef mutable struct StormChaser{NF} <: SpeedyWeather.AbstractCallback
        timestep_counter::Int = 0
        maximum_surface_wind_speed::Vector{NF} = Float64[]
    end
    
    # Generator function
    StormChaser(SG::SpectralGrid) = StormChaser{SG.NF}()

    function SpeedyWeather.initialize!(
        callback::StormChaser,
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        model::ModelSetup,
    )
        # allocate recorder: number of time steps (incl initial conditions) in simulation  
        callback.maximum_surface_wind_speed = zeros(progn.clock.n_timesteps + 1)
        
        # where surface (=lowermost model layer) u,v on the grid are stored
        (;u_grid, v_grid) = diagn.layers[diagn.nlev].grid_variables
        
        # maximum wind speed of initial conditions
        callback.maximum_surface_wind_speed[1] = max_2norm(u_grid,v_grid)
        
        # (re)set counter to 1
        callback.timestep_counter = 1
    end

    function max_2norm(u::AbstractArray{T},v::AbstractArray{T}) where T
        max_norm = zero(T)      # = u² + v²
        for ij in eachindex(u, v)
            # find largest wind speed squared
            max_norm = max(max_norm, u[ij]^2 + v[ij]^2)
        end
        return sqrt(max_norm)   # take sqrt only once
    end

    function SpeedyWeather.callback!(
        callback::StormChaser,
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        model::ModelSetup,
    )
    
        # increase counter
        callback.timestep_counter += 1  
        i = callback.timestep_counter
    
        # where surface (=lowermost model layer) u,v on the grid are stored
        (;u_grid, v_grid) = diagn.layers[diagn.nlev].grid_variables
    
        # maximum wind speed at current time step
        callback.maximum_surface_wind_speed[i] = max_2norm(u_grid,v_grid)
    end

    SpeedyWeather.finish!(::StormChaser,args...) = nothing

    spectral_grid = SpectralGrid()
    # callbacks = [NoCallback()]    # doesn't work at the moment
    model = PrimitiveWetModel(;spectral_grid)
    
    storm_chaser = StormChaser(spectral_grid)
    append!(model.callbacks, storm_chaser)
    append!(model.callbacks, NoCallback())  # add dummy too 

    simulation = initialize!(model)
    run!(simulation)

    # maximum wind speed should always be non-negative
    @test all(model.callbacks[1].maximum_surface_wind_speed .>= 0)

    # highest wind speed across all time steps should be positive
    @test maximum(model.callbacks[1].maximum_surface_wind_speed) > 0
end