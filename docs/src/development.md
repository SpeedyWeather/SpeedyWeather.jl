# Development notes

To run tests, from the path of your local clone of the repository do:
```
julia --project=. -e 'import Pkg; Pkg.test()'
```

To install dependencies:

```
julia --project -e 'import Pkg; Pkg.instantiate()`
```

then:
```
julia --project=.
```

To generate docs:

```
julia --project=docs -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
