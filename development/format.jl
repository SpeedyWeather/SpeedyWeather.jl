import Pkg

Pkg.instantiate()

using JuliaFormatter

println("Running JuliaFormatter...")
format(".", verbose = false, overwrite = true)
