# Methodology

## 1. Keplerian propagation

The two-body model uses

`n = sqrt(μ/a^3)`

and

`M(t) = M₀ + nt`.

The eccentric anomaly is obtained from Kepler's equation,

`M = E - e sin(E)`,

using Newton–Raphson iteration.

The solved eccentric anomaly is converted to a PQW position and velocity,
then transformed into ECI coordinates using the classical
`R3(-Ω) R1(-i) R3(-ω)` sequence.

## 2. Earth-fixed coordinates

The implementation applies the simplified constant-rate Earth rotation model

`GAST = ωE t`

and rotates the ECI position into ECEF coordinates. Geocentric latitude and
longitude are then computed from the spherical-Earth position vector.

## 3. Secular J₂ model

The implemented secular rates are

`Ω̇ = -(3/2)n J₂ (Re/a)^2 cos(i)/(1-e²)^2`

`ω̇ = (3/4)n J₂ (Re/a)^2 (5cos²(i)-1)/(1-e²)^2`

`Ṁ = n + (3/4)n J₂ (Re/a)^2 (3cos²(i)-1)/(1-e²)^(3/2)`.

The model keeps `a`, `e`, and `i` constant and evolves Ω, ω and M linearly.

## 4. Repeat-orbit design

The repeat condition is formulated as

`(Ṁ + ω̇)/(ωE - Ω̇) = p/q`.

The design problem is therefore reduced to a scalar root-finding problem in
semi-major axis `a`. Bisection is used after verifying that the requested root
is bracketed.

The API explicitly names `p` as orbital revolutions and `q` as Earth rotations,
avoiding ambiguity in the original α/β notation.

## 5. Closure analysis

For a selected repeat period, a reference ECEF trajectory is sampled over one
repeat cycle. The same sample epochs are shifted by integer multiples of the
repeat period and compared against the reference.

For each cycle the framework reports:

- mean position error,
- RMS position error,
- maximum position error.

This provides a numerical repeatability metric in addition to visual ground-track
inspection.

## 6. Planned independent numerical validation

The remaining extension is a numerical orbit integrator using an independent
ODE formulation and force model. The intended comparison should include:

- ECI position residuals,
- ECI velocity residuals,
- orbital-element differences,
- ground-track differences,
- growth of residuals with propagation time.

That comparison should be added before claiming full numerical-integrator
validation in a CV or project report.
