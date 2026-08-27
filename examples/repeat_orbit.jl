include(joinpath(@__DIR__, "..", "src", "SatelliteOrbitPropagation.jl"))
using .SatelliteOrbitPropagation

# Example: 15 orbital revolutions per 1 Earth rotation.
p = 15
q = 1
inclination = deg2rad(98.0)

a = design_repeat_orbit(
    p,
    q,
    inclination;
    e = 0.0,
)

println("Designed semi-major axis = ", a, " km")
println("J2-aware repeat ratio    = ",
        repeat_ratio(a, 0.0, inclination))

sat = Satellite(
    "15x1 J2 Repeat Orbit",
    a,
    0.0,
    inclination,
    0.0,
    0.0,
    0.0,
)

T_orbit = orbital_period(a)
T_repeat = repeat_period(sat, p, q)

println("Orbital period = ", T_orbit / 60, " min")
println("Repeat period  = ", T_repeat / 3600, " hr")

result = propagate_j2(
    sat;
    dt = 30.0,
    revolutions = 30,
    repeat_revolutions = p,
)

summary(result)
display(plot_groundtrack(result; title = "15/1 J2 Repeat Ground Track"))

closure = repeat_cycle_closure(
    sat,
    T_repeat;
    n_cycles = 5,
    n_samples = 500,
)

println()
println("Repeat-cycle closure statistics")
for k in eachindex(closure.rms_error)
    println(
        "Cycle ", k,
        ": RMS = ", closure.rms_error[k],
        " km, max = ", closure.maximum_error[k], " km"
    )
end

rate_check = element_rate_check(sat, T_repeat)
println()
println("Secular-rate residuals at one repeat period:")
println(rate_check.residual)
