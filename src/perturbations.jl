function compute_j2_secular_rates(
    a::Real,
    e::Real,
    i::Real;
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
)
    a > 0 || throw(ArgumentError("Semi-major axis must be positive."))
    0 <= e < 1 || throw(ArgumentError("Eccentricity must satisfy 0 ≤ e < 1."))

    n = sqrt(μ / a^3)
    factor = J2 * (Re / a)^2
    β = 1 - e^2
    ci = cos(i)

    Ωdot = -(3 / 2) * n * factor * ci / β^2
    ωdot =  (3 / 4) * n * factor * (5ci^2 - 1) / β^2
    Mdot = n + (3 / 4) * n * factor * (3ci^2 - 1) / β^(3 / 2)

    return Ωdot, ωdot, Mdot
end

function propagate_j2_elements(
    sat::Satellite,
    t::Real;
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
)
    Ωdot, ωdot, Mdot = compute_j2_secular_rates(
        sat.a, sat.e, sat.i;
        μ = μ, Re = Re, J2 = J2,
    )

    return Satellite(
        sat.name,
        sat.a,
        sat.e,
        sat.i,
        sat.Ω + Ωdot * t,
        sat.ω + ωdot * t,
        sat.M0 + Mdot * t,
    )
end
