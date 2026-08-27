using Plots

"""
    plot_orbital_elements(result)

Plots the evolution of all six classical orbital elements.

Input:
    result : PropagationResult

Returns:
    A vector containing the six plots.
"""
function plot_orbital_elements(result::PropagationResult)

    t = result.time ./ 3600      # hours

    e = result.elements

    p1 = plot(
        t,
        e.a,
        xlabel="Time (hr)",
        ylabel="a (km)",
        title="Semi-major Axis",
        legend=false
    )

    p2 = plot(
        t,
        e.e,
        xlabel="Time (hr)",
        ylabel="e",
        title="Eccentricity",
        legend=false
    )

    p3 = plot(
        t,
        rad2deg.(e.i),
        xlabel="Time (hr)",
        ylabel="i (deg)",
        title="Inclination",
        legend=false
    )

    p4 = plot(
        t,
        rad2deg.(e.Ω),
        xlabel="Time (hr)",
        ylabel="Ω (deg)",
        title="RAAN",
        legend=false
    )

    p5 = plot(
        t,
        rad2deg.(e.ω),
        xlabel="Time (hr)",
        ylabel="ω (deg)",
        title="Argument of Perigee",
        legend=false
    )

    p6 = plot(
        t,
        rad2deg.(e.M),
        xlabel="Time (hr)",
        ylabel="M (deg)",
        title="Mean Anomaly",
        legend=false
    )

    return [p1,p2,p3,p4,p5,p6]

end