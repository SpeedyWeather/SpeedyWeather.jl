abstract type AbstractTimeStepper <: AbstractModelComponent end
abstract type AbstractImplicit <: AbstractModelComponent end

"""Placeholder for the implicit solver of a time stepper, replaced by the model's default
(`default_implicit`) when the model is constructed, e.g. `nothing` for the `BarotropicModel`,
`ImplicitShallowWater` for the `ShallowWaterModel` and `ImplicitPrimitiveEquation` for the
primitive equation models."""
struct DefaultImplicit <: AbstractImplicit end
abstract type AbstractDynamicalCoreComponent <: AbstractModelComponent end

# DUMMY TYPES
# these are used for dispatching when get_..._step but the place
# where isn't a model.component per se, so we define dummy types
# that can be used instead to distinguish and let the time stepper
# decide which step to get

struct DummyTimeStepper <: AbstractTimeStepper end
struct DummyParameterization <: AbstractParameterization end
struct DynamicalCore <: AbstractDynamicalCoreComponent end
struct LinearDynamicalCore <: AbstractDynamicalCoreComponent end
struct BernoulliPotential <: AbstractDynamicalCoreComponent end
struct ContinuityEquation <: AbstractDynamicalCoreComponent end
struct ResetTendencies <: AbstractDynamicalCoreComponent end
