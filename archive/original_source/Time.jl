"""
function mean_motion(a, GM)

Computes the mean motion of a Keplerian orbit.

Inputs:
a  : Semi-major axis [km]
GM : Gravitational parameter [km^3/s^2]

Returns:
n : Mean motion [rad/s]


function orbital_period(a, GM)

Computes the orbital period of a Keplerian orbit.

Inputs:
a  : Semi-major axis [km]
GM : Gravitational parameter [km^3/s^2]

Returns:
T : Orbital period [s]


function mean_anomaly(t, n, M0 = 0.0)

Computes the mean anomaly at time t.

Inputs:
t  : Time since epoch [s]
n  : Mean motion [rad/s]
M0 : Mean anomaly at epoch [rad]

Returns:
M : Mean anomaly [rad]


function gast(t, ωE)

Computes the Greenwich Apparent Sidereal Time (GAST)
assuming a constant Earth rotation rate.

Inputs:
t  : Time since epoch [s]
ωE : Earth's angular velocity [rad/s]

Returns:
GAST : Greenwich Apparent Sidereal Time [rad]

"""
function mean_motion(a, GM)

    n = sqrt(GM / a^3)

    return n

end

function orbital_period(a, GM)

    n = mean_motion(a, GM)

    T = 2π / n

    return T

end

function mean_anomaly(t, n, M0 = 0.0)

    M = M0 + n * t

    return M

end


function gast(t, ωE)

    GAST = ωE * t

    return GAST

end