using ParallelTestRunner
using SpeedyWeather

# code that's run for every worker before running the tests
const init_code = quote
    using SpeedyWeather
end

# Start with autodiscovered tests
testsuite = find_tests(@__DIR__)

# Parse arguments
args = parse_args(ARGS)

# We don't run the following tests using `Pkg.test`
delete!(testsuite, "parameters")
delete!(testsuite, "prognostic_variables_stability_test")
delete!(testsuite, "type_stability_test")
for key in keys(testsuite)
    if startswith(key, r"(GPU|differentiability|reactant)/") && key != "GPU/kernelabstractions"
        delete!(testsuite, key)
    end
end

# long integration tests should always be run with `-O2` compiler flag, as set in workflow yml,
# return nothing uses default worker, but all other tests should use `-O0` as compile time heavy
test_worker(name) = startswith(name, "long") ? nothing : addworker(; exeflags = ["-O0"])

# run tests in parallel
runtests(SpeedyWeather, args; test_worker, testsuite, init_code)
