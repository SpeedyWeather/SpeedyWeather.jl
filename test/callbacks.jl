@testset "CallbackDict generator" begin
    @test CallbackDict() isa SpeedyWeather.CALLBACK_DICT
    @test CallbackDict(NoCallback()) isa SpeedyWeather.CALLBACK_DICT
    @test CallbackDict(GlobalSurfaceTemperatureCallback()) isa SpeedyWeather.CALLBACK_DICT
    @test CallbackDict(NoCallback(), NoCallback()) isa SpeedyWeather.CALLBACK_DICT
    @test CallbackDict(NoCallback(), GlobalSurfaceTemperatureCallback()) isa SpeedyWeather.CALLBACK_DICT
    @test CallbackDict(:a => NoCallback()) isa SpeedyWeather.CALLBACK_DICT

    d = CallbackDict()
    add!(d, NoCallback())
    add!(d, NoCallback(), NoCallback())
    add!(d, :my_callback, NoCallback())
    add!(d, "my_callback2", NoCallback())
    delete!(d, :my_callback)

    add!(d, :another_callback => NoCallback())
    add!(d, :another_callback1 => NoCallback(), :another_callback2 => NoCallback())
    @test d isa SpeedyWeather.CALLBACK_DICT
    @test length(d) == 7
end

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
        model::AbstractModel,
    )
        # allocate recorder: number of time steps (incl initial conditions) in simulation  
        callback.maximum_surface_wind_speed = zeros(progn.clock.n_timesteps + 1)
        
        # where surface (=lowermost model layer) u, v on the grid are stored
        u_grid = diagn.grid.u_grid[:, diagn.nlayers]
        v_grid = diagn.grid.u_grid[:, diagn.nlayers]

        # maximum wind speed of initial conditions
        callback.maximum_surface_wind_speed[1] = max_2norm(u_grid, v_grid)
        
        # (re)set counter to 1
        callback.timestep_counter = 1
    end

    function max_2norm(u::AbstractArray{T}, v::AbstractArray{T}) where T
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
        model::AbstractModel,
    )
    
        # increase counter
        callback.timestep_counter += 1  
        i = callback.timestep_counter
    
        # where surface (=lowermost model layer) u, v on the grid are stored
        u_grid = diagn.grid.u_grid[:,diagn.nlayers]
        v_grid = diagn.grid.u_grid[:,diagn.nlayers]
        
        # maximum wind speed at current time step
        callback.maximum_surface_wind_speed[i] = max_2norm(u_grid, v_grid)
    end

    SpeedyWeather.finalize!(::StormChaser, args...) = nothing

    spectral_grid = SpectralGrid()
    callbacks = CallbackDict(NoCallback())
    model = PrimitiveWetModel(spectral_grid; callbacks)
    
    storm_chaser = StormChaser(spectral_grid)
    key = :storm_chaser
    add!(model.callbacks, key => storm_chaser)      # with :storm_chaser key
    add!(model.callbacks, NoCallback())             # add dummy too 
    add!(model, NoCallback())                       # add dummy with ::AbstractModel interface

    simulation = initialize!(model)
    run!(simulation, period=Day(1))

    # maximum wind speed should always be non-negative
    @test all(model.callbacks[key].maximum_surface_wind_speed .>= 0)

    # highest wind speed across all time steps should be positive
    @test maximum(model.callbacks[key].maximum_surface_wind_speed) > 0
end