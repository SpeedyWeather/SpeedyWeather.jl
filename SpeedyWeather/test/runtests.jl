using ParallelTestRunner
using SpeedyWeather

# code that's run for every worker before running the tests
const init_code = quote
    using SpeedyWeather
end

# test suites, manual or automatic file discovery
testsuite_GPU = Dict(
    "kernelabstractions" => quote
        include("GPU/kernelabstractions.jl")
    end
)

testsuite_dynamics = find_tests(joinpath(pwd(), "dynamics"))
testsuite_long_integrations = find_tests(joinpath(pwd(), "long_integrations.jl"))
testsuite_parameterizations = find_tests(joinpath(pwd(), "parameterizations"))
testsuite_output = find_tests(joinpath(pwd(), "output"))
testsuite_variables = find_tests(joinpath(pwd(), "variables"))

# merge all testsuites
testsuite = merge(
    testsuite_GPU,
    testsuite_dynamics,
    testsuite_long_integrations,
    testsuite_parameterizations,
    testsuite_variables,
    testsuite_output
)

# long integration tests should always be run with `-O2` compiler flag, as set in workflow yml,
# return nothing uses default worker, but all other tests should use `-O0` as compile time heavy
# 
test_worker(name) = contains(name, "long") ? nothing : addworker(; exeflags = ["-O0"])

# run tests in parallel
runtests(SpeedyWeather, ARGS; test_worker, testsuite, init_code)
