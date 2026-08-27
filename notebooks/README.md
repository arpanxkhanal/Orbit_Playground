# Notebooks

The original exploratory notebooks supplied during development are kept under
`archive/notebooks/`.

The maintained project workflow is intentionally executable from the scripts in
`examples/`, so the repository does not depend on a notebook runtime for its
core functionality.

A future polished notebook can import the package with:

```julia
include("../src/SatelliteOrbitPropagation.jl")
using .SatelliteOrbitPropagation
```
