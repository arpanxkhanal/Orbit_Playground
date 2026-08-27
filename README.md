# Satellite Orbit Propagation & Repeat-Orbit Design Framework

A modular Julia framework for analytical satellite-orbit propagation, Earth-fixed
ground-track analysis, secular J₂ perturbations, and repeat-orbit design.

The project was developed as an independent computational project in satellite
orbit mechanics and geoinformatics, with an emphasis on reusable numerical
components rather than a single script.

## What is implemented

- Two-body Keplerian orbit propagation.
- Newton–Raphson solution of Kepler's equation with convergence control.
- State construction in the perifocal/PQW frame.
- PQW → ECI and ECI → ECEF coordinate transformations.
- Geocentric latitude/longitude and ground-track generation.
- Analytical secular J₂ evolution of Ω, ω, and mean anomaly.
- Time histories of all six classical orbital elements.
- Repeat-orbit condition including J₂ secular rates.
- Bisection-based semi-major-axis design.
- Multi-revolution repeat-orbit propagation.
- Quantitative ECEF repeat-cycle closure analysis using RMS and maximum errors.
- Plotting utilities for 3-D orbits, ground tracks, and orbital-element histories.
- Unit tests for the main analytical components.

## Scope and limitations

The analytical J₂ model used here is a **secular/mean-element approximation**:
`a`, `e`, and `i` remain constant while Ω, ω, and M evolve linearly according
to the implemented secular rates. It is not a high-fidelity force model.

The Earth-fixed transformation uses a constant Earth rotation rate and a simplified
GAST model. Latitude/longitude are geocentric and use a spherical Earth.

### Not yet implemented

The remaining planned validation item is comparison against an independent
numerical orbit integrator, including state residuals and long-duration orbital
behaviour. A clean extension point is reserved for this in `src/validation.jl`.

## Repeat-orbit convention

The public API avoids the ambiguous α/β notation used in the original scripts.
A repeat ratio is written explicitly as:

`orbits_per_repeat / earth_rotations_per_repeat`

For example, a conventional **15/1 repeat orbit** means approximately 15 orbital
revolutions per one Earth rotation over the repeat cycle.

The J₂-aware condition implemented by `repeat_ratio()` is

`(Ṁ + ω̇) / (ωE - Ω̇) = p / q`

where `p` is the number of orbital revolutions and `q` is the number of Earth
rotations in the repeat cycle.

## Quick start

Install Julia 1.9+ and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then:

```julia
include("src/SatelliteOrbitPropagation.jl")
using .SatelliteOrbitPropagation

sat = Satellite(
    "Demo",
    7000.0,
    0.01,
    deg2rad(98.0),
    0.0,
    0.0,
    0.0,
)

result = propagate_kepler(
    sat;
    revolutions = 3,
    dt = 20.0,
)

plot_groundtrack(result)
```

For J₂ propagation:

```julia
result_j2 = propagate_j2(
    sat;
    revolutions = 15,
    dt = 30.0,
)
```

For repeat-orbit design:

```julia
a = design_repeat_orbit(
    15,
    1,
    deg2rad(98.0);
    e = 0.0,
)

repeat = Satellite(
    "15/1 repeat candidate",
    a,
    0.0,
    deg2rad(98.0),
    0.0,
    0.0,
    0.0,
)

closure = repeat_cycle_closure(
    repeat,
    repeat_period(repeat, 15, 1);
    n_cycles = 5,
    n_samples = 500,
)
```

See `examples/repeat_orbit.jl` for a complete workflow.

## Repository layout

```text
satellite-orbit-propagation/
├── Project.toml
├── README.md
├── LICENSE
├── src/
│   ├── SatelliteOrbitPropagation.jl
│   ├── constants.jl
│   ├── types.jl
│   ├── kepler.jl
│   ├── time.jl
│   ├── rotations.jl
│   ├── coordinates.jl
│   ├── perturbations.jl
│   ├── propagation.jl
│   ├── repeat_orbit.jl
│   ├── analysis.jl
│   ├── plotting.jl
│   └── validation.jl
├── examples/
│   ├── propagation_demo.jl
│   └── repeat_orbit.jl
├── test/
│   └── runtests.jl
├── notebooks/
│   └── README.md
├── archive/
│   ├── original_source/
│   └── notebooks/
└── docs/
    └── methodology.md
```

## Validation philosophy

The test suite checks mathematical identities and internal consistency:
Kepler residuals, orbital-period/mean-motion consistency, orthogonality of
rotation matrices, state construction, J₂ secular rates, repeat-condition
residuals, bisection convergence, and closure-analysis dimensions.

The numerical-integrator comparison is deliberately not represented as complete
until an independent numerical propagator is added.

## Placement-project framing

This repository is intended to present the work as a small research/computational
framework rather than as a collection of exploratory Julia files. The original
scripts and notebooks are retained under `archive/` for traceability, while the
maintained implementation lives under `src/`.
