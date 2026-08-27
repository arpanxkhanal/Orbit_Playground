function solve_kepler(
    e::Real,
    M::Real;
    tol::Real = 1e-12,
    max_iter::Int = 50,
)
    0 <= e < 1 || throw(ArgumentError("Eccentricity must satisfy 0 ≤ e < 1."))
    tol > 0 || throw(ArgumentError("Tolerance must be positive."))
    max_iter > 0 || throw(ArgumentError("max_iter must be positive."))

    # Second-order initial estimate followed by Newton–Raphson iterations.
    E = M + e * sin(M) + 0.5 * e^2 * sin(2M)

    for _ in 1:max_iter
        f = E - e * sin(E) - M
        dfdE = 1 - e * cos(E)
        ΔE = -f / dfdE
        E += ΔE

        if abs(ΔE) < tol
            return E
        end
    end

    throw(ErrorException(
        "Kepler's equation did not converge after $max_iter iterations."
    ))
end
