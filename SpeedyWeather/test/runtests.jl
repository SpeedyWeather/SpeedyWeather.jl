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
testsuite_parameterizations = find_tests(joinpath(pwd(), "parameterizations"))
testsuite_output = find_tests(joinpath(pwd(), "output"))

# merge all testsuites
testsuite = merge(
    testsuite_GPU,
    testsuite_dynamics,
    testsuite_parameterizations,
    testsuite_output
)

# run tests in parallel
runtests(SpeedyWeather, ARGS; testsuite, init_code)
