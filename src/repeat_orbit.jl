function repeat_ratio(
    a::Real,
    e::Real,
    i::Real;
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
    ωE::Real = EARTH_ROTATION_RATE,
)
    Ωdot, ωdot, Mdot = compute_j2_secular_rates(
        a, e, i;
        μ = μ, Re = Re, J2 = J2,
    )

    return (Mdot + ωdot) / (ωE - Ωdot)
end

function repeat_orbit_error(
    a::Real,
    e::Real,
    i::Real,
    p::Integer,
    q::Integer;
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
    ωE::Real = EARTH_ROTATION_RATE,
)
    p > 0 || throw(ArgumentError("p must be positive."))
    q > 0 || throw(ArgumentError("q must be positive."))

    lhs = repeat_ratio(
        a, e, i;
        μ = μ, Re = Re, J2 = J2, ωE = ωE,
    )
    rhs = p / q

    return lhs - rhs
end

function design_repeat_orbit(
    p::Integer,
    q::Integer,
    i::Real;
    e::Real = 0.0,
    a_low::Real = 6500.0,
    a_high::Real = 8000.0,
    tol::Real = 1e-10,
    max_iter::Int = 100,
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
    ωE::Real = EARTH_ROTATION_RATE,
)
    p > 0 || throw(ArgumentError("p must be positive."))
    q > 0 || throw(ArgumentError("q must be positive."))
    a_low < a_high || throw(ArgumentError("a_low must be less than a_high."))
    tol > 0 || throw(ArgumentError("tol must be positive."))

    f(a) = repeat_orbit_error(
        a, e, i, p, q;
        μ = μ, Re = Re, J2 = J2, ωE = ωE,
    )

    f_low = f(a_low)
    f_high = f(a_high)

    f_low * f_high <= 0 || throw(ArgumentError(
        "Repeat-orbit root is not bracketed. " *
        "Try different semi-major-axis bounds."
    ))

    for _ in 1:max_iter
        a_mid = (a_low + a_high) / 2
        f_mid = f(a_mid)

        if abs(f_mid) < tol
            return a_mid
        end

        if f_low * f_mid <= 0
            a_high = a_mid
            f_high = f_mid
        else
            a_low = a_mid
            f_low = f_mid
        end
    end

    throw(ErrorException("Bisection did not converge within max_iter."))
end

function repeat_period(
    sat::Satellite,
    p::Integer,
    q::Integer;
    μ::Real = EARTH_MU,
)
    a = sat.a
    return p * orbital_period(a, μ)
end
