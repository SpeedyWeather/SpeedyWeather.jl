struct DummyArchitecture end
const DEFAULT_ARCHITECTURE = DummyArchitecture
array_type(::Type{DummyArchitecture}) = Array
array_type(arch::DummyArchitecture) = array_type(typeof(arch))