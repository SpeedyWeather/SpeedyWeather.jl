abstract type AbstractTimeStepper <: AbstractModelComponent end
abstract type AbstractImplicit <: AbstractModelComponent end
abstract type AbstractDynamicalCoreComponent <: AbstractModelComponent end

struct DummyTimeStepper <: AbstractTimeStepper end

# define dummy structs for components within the dynamical core
# to enable dispatch with get_step for them
struct DynamicalCore <: AbstractDynamicalCoreComponent end
struct BernoulliPotential <: AbstractDynamicalCoreComponent end
struct ContinuityEquation <: AbstractDynamicalCoreComponent end
struct ResetTendencies <: AbstractDynamicalCoreComponent end
