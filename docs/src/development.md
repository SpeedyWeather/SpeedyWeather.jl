# Development notes

To run tests, from the path of your local clone of the repository do:

```
julia --project=. -e 'import Pkg; Pkg.test()'
```

To install dependencies:

```
julia --project -e 'import Pkg; Pkg.instantiate()`
```

then opening:

```
julia --project
```

you are able to `using SpeedyWeather`.

To generate the docs:

```
julia --project=docs -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

If the docs are generated successfully, you view them by opening `docs/build/index.html` in your favorite browser.
