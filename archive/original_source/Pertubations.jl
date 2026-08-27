"""
Perturbations.jl

This file contains functions related to analytical orbit perturbations.

Currently Implemented:
- J₂ Secular Perturbation
"""

"""
compute_j2_secular_rates(
    a,
    e,
    i;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3
)

Computes the secular drift rates of the Keplerian elements
due to the Earth's J₂ perturbation.

Inputs:
a : Semi-major axis [km]
e : Eccentricity
i : Inclination [rad]

Keyword Arguments:
μ  : Earth's gravitational parameter [km³/s²]
Re : Earth's equatorial radius [km]
J2 : Second zonal harmonic coefficient

Returns:
Ω̇ : RAAN rate [rad/s]
ω̇ : Argument of Perigee rate [rad/s]
Ṁ : Mean Anomaly rate [rad/s]
"""
function compute_j2_secular_rates(
    a,
    e,
    i;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3
)
    # Mean motion
    n = sqrt(μ / a^3)

    # Common J₂ factor
    factor = J2 * (Re / a)^2

    # Eccentricity term
    β = 1.0 - e^2

    # Cosine of inclination
    ci = cos(i)

    Ωdot = -(3/2) * n * factor * ci / β^2

    ωdot = (3/4) * n * factor * (5 * ci^2 - 1) / β^2

    Mdot = n + (3/4) * n * factor * (3 * ci^2 - 1) / β^(3/2)

    return Ωdot, ωdot, Mdot
end




"""
propagate_j2_elements(
    sat,
    t;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3
)

Propagates the Keplerian orbital elements using the
analytical secular J₂ perturbation model.

Only the following elements are updated:
- Right Ascension of the Ascending Node (Ω)
- Argument of Perigee (ω)
- Mean Anomaly (M₀)

The remaining elements (a, e, i) are assumed constant.

Inputs:
sat : Satellite structure
t   : Time since epoch [s]

Keyword Arguments:
μ  : Earth's gravitational parameter [km³/s²]
Re : Earth's equatorial radius [km]
J2 : Second zonal harmonic coefficient

Returns:
Satellite : Updated satellite with propagated orbital elements.
"""

function propagate_j2_elements(
    sat::Satellite,
    t;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3
)

    # Compute secular drift rates
    Ωdot, ωdot, Mdot = compute_j2_secular_rates(
        sat.a,
        sat.e,
        sat.i;
        μ,
        Re,
        J2
    )

    # Update Keplerian elements
    Ω  = sat.Ω  + Ωdot * t
    ω  = sat.ω  + ωdot * t
    M0 = sat.M0 + Mdot * t

    # Return updated satellite
    N = length(t_list)
    return Satellite(
        sat.name,
        sat.a,
        sat.e,
        sat.i,
        Ω,
        ω,
        M0
    )

end