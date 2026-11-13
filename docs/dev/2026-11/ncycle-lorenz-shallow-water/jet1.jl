include("gw.jl")
WDIFF[] = get(ENV, "WDIFF", "0") == "1"
trunc = parse(Int, ARGS[1]); dt = parse(Int, ARGS[2]); θs = ARGS[3]
sg = SpectralGrid(truncation=trunc, nlayers=1)
time_stepping = θs == "LF" ? Leapfrog(sg; Δt_at_T32=Minute(dt), adjust_with_output=false) : NCycleLorenz(sg; Δt_at_T32=Minute(dt), steps=3, adjust_with_output=false)
implicit = θs in ("LF","old") ? ImplicitShallowWater(sg) : θs == "none" ? nothing : ImplicitSWθ{Float32}(Float32.(parse.(Float64, split(θs, ","))))
model = ShallowWaterModel(sg; time_stepping, implicit, initial_conditions=ZonalJet(sg))
sim = initialize!(model)
ok = try run!(sim, period=Day(20)); !model.feedback.nans_detected catch; false end
@printf("T%d Δt=%d %-22s stable=%s\n", trunc, dt, θs, ok)
