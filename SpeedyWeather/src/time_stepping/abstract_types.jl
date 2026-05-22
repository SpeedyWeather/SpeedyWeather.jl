abstract type AbstractTimeStepper <: AbstractModelComponent end
abstract type AbstractImplicit <: AbstractModelComponent end
abstract type AbstractDynamicalCoreComponent <: AbstractModelComponent end

# DUMMY TYPES
# these are used for dispatching when get_..._step but the place
# where isn't a model.component per se, so we define dummy types
# that can be used instead to distinguish and let the time stepper
# decide which step to get

struct DummyTimeStepper <: AbstractTimeStepper end
struct DummyParameterization <: AbstractParameterization end
struct DynamicalCore <: AbstractDynamicalCoreComponent end
struct BernoulliPotential <: AbstractDynamicalCoreComponent end
struct ContinuityEquation <: AbstractDynamicalCoreComponent end
struct ResetTendencies <: AbstractDynamicalCoreComponent end