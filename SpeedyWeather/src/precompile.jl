# Precompile the setup and the first time steps of the default model, with and without
# NetCDF output, so that `initialize!(model)` and `run!(simulation)` are fast in a new session.
# Compilation is specific to the number format: only NF = Float32 (the default) is covered.
# Truncation and nlayers can be changed without recompiling (most of) the model.
# Disable this workload (e.g. for faster precompilation during development) with
#   using SpeedyWeather, Preferences
#   Preferences.set_preferences!(SpeedyWeather, "precompile_workload" => false; force = true)
@setup_workload begin
    path = mktempdir()
    @compile_workload begin
        # the default model reads boundary conditions from assets that may need to be downloaded
        # first, if that fails (e.g. no internet) skip the workload rather than failing precompilation
        try
            spectral_grid = SpectralGrid()
            output = NetCDFOutput(spectral_grid; path)
            model = PrimitiveWetModel(spectral_grid; output)
            simulation = initialize!(model)
            run!(simulation, steps = 2, output = true)
        catch e
            @warn "SpeedyWeather precompile workload failed, compilation deferred to first use." exception = e
        end
    end
    rm(path, recursive = true, force = true)
end
