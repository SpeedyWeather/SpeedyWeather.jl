# The drag's neutral transfer coefficient has a finite derivative w.r.t. the height of the lowermost
# layer also for an absent surface type (zero roughness length): log(z/0) = Inf gave a NaN derivative
# (0 · ∞) in reverse mode with Reactant, which made the gradient of the PrimitiveWetModel non-finite

coefficient_sum(z, z₀) = sum(SpeedyWeather.neutral_transfer_coefficient.(z, z₀, z₀, 0.4f0))

@testset "Neutral transfer coefficient derivative with zero roughness (Reactant)" begin
    z, dz = Reactant.to_rarray([30.0f0, 30.0f0]), Reactant.to_rarray([0.0f0, 0.0f0])
    z₀, dz₀ = Reactant.to_rarray([0.0f0, 1.0f-3]), Reactant.to_rarray([0.0f0, 0.0f0])   # absent and present surface
    @jit Enzyme.autodiff(Reverse, coefficient_sum, Active, Duplicated(z, dz), Duplicated(z₀, dz₀))
    @test all(isfinite, Array(dz))
    @test all(isfinite, Array(dz₀))
end
