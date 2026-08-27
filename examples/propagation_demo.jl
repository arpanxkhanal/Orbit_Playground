include(joinpath(@__DIR__, "..", "src", "SatelliteOrbitPropagation.jl"))
using .SatelliteOrbitPropagation

sat = Satellite(
    "DemoSat",
    7000.0,
    0.01,
    deg2rad(98.0),
    0.0,
    0.0,
    0.0,
)

println("Running Keplerian propagation...")
kepler_result = propagate_kepler(
    sat;
    dt = 20.0,
    revolutions = 3,
)
summary(kepler_result)

display(plot_orbit(kepler_result.positions_eci))
display(plot_groundtrack(kepler_result))

println("Running analytical secular J2 propagation...")
j2_result = propagate_j2(
    sat;
    dt = 20.0,
    revolutions = 3,
)
summary(j2_result)

display(plot_groundtrack(j2_result))
display(plot_orbital_elements(j2_result)[4])
