"""
repeat_condition(
    a,
    e,
    i,
    α,
    β;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3,
    ωE = 7.292115e-5
)

Evaluates the repeat orbit condition

    (n + ω̇)/(ωE − Ω̇) = β/α

Returns

lhs : Computed repeat ratio
rhs : Desired repeat ratio
error : lhs - rhs
"""
function repeat_condition(
    a,
    e,
    i,
    α,
    β;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3,
    ωE = 7.292115e-5
)

    # Mean motion
    n = sqrt(μ / a^3)

    # J2 secular rates
    Ωdot, ωdot, Mdot = compute_j2_secular_rates(
        a,
        e,
        i;
        μ = μ,
        Re = Re,
        J2 = J2
    )

    lhs = (Mdot + ωdot) / (ωE - Ωdot)

    rhs = β / α

    return lhs, rhs, lhs - rhs

end

"""
repeat_orbit_error(
    a,
    e,
    i,
    α,
    β;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3,
    ωE = 7.292115e-5
)

Returns the repeat orbit error.

A perfect repeat orbit satisfies

error = 0
"""
function repeat_orbit_error(
    a,
    e,
    i,
    α,
    β;
    μ = 398600.4418,
    Re = 6378.137,
    J2 = 1.08263e-3,
    ωE = 7.292115e-5
)

    lhs, rhs, err = repeat_condition(
        a,
        e,
        i,
        α,
        β;
        μ = μ,
        Re = Re,
        J2 = J2,
        ωE = ωE
    )

    return err

end

"""
design_repeat_orbit(
    α,
    β,
    i;
    e = 0.0,
    a_low = 6500.0,
    a_high = 8000.0,
    tol = 1e-10,
    max_iter = 100
)

Designs a repeat orbit by solving for the semi-major axis.
"""
function design_repeat_orbit(
    α,
    β,
    i;
    e = 0.0,
    a_low = 400.0,
    a_high = 800.0,
    tol = 1e-10,
    max_iter = 100
)

    err_low = repeat_orbit_error(a_low, e, i, α, β)
    err_high = repeat_orbit_error(a_high, e, i, α, β)

    if err_low * err_high > 0
        error("Root is not bracketed. Choose different bounds.")
    end

    for k in 1:max_iter

        a_mid = (a_low + a_high) / 2

        err_mid = repeat_orbit_error(a_mid, e, i, α, β)

        if abs(err_mid) < tol
            return a_mid
        end

        if err_low * err_mid < 0
            a_high = a_mid
            err_high = err_mid
        else
            a_low = a_mid
            err_low = err_mid
        end

    end

    error("Bisection did not converge.")

end