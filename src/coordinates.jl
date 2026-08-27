function kepler_to_pqw(sat::Satellite, E::Real, μ::Real = EARTH_MU)
    a = sat.a
    e = sat.e

    a > 0 || throw(ArgumentError("Semi-major axis must be positive."))
    0 <= e < 1 || throw(ArgumentError("Eccentricity must satisfy 0 ≤ e < 1."))

    x = a * (cos(E) - e)
    y = a * sqrt(1 - e^2) * sin(E)
    r_pqw = [x, y, 0.0]

    r = norm(r_pqw)
    factor = sqrt(μ * a) / r

    v_pqw = [
        -factor * sin(E),
         factor * sqrt(1 - e^2) * cos(E),
         0.0,
    ]

    return r_pqw, v_pqw
end

function pqw_to_eci(
    sat::Satellite,
    r_pqw::AbstractVector,
    v_pqw::AbstractVector,
)
    length(r_pqw) == 3 || throw(ArgumentError("r_pqw must have length 3."))
    length(v_pqw) == 3 || throw(ArgumentError("v_pqw must have length 3."))

    Q = R3(-sat.Ω) * R1(-sat.i) * R3(-sat.ω)
    return Q * r_pqw, Q * v_pqw
end

function eci_to_ecef(
    r_eci::AbstractVector,
    θ::Real,
)
    length(r_eci) == 3 || throw(ArgumentError("r_eci must have length 3."))
    return R3(θ) * r_eci
end

function ecef_to_latlon(r_ecef::AbstractVector)
    length(r_ecef) == 3 || throw(ArgumentError("r_ecef must have length 3."))

    r = norm(r_ecef)
    r > 0 || throw(ArgumentError("Position vector must be non-zero."))

    x, y, z = r_ecef
    lat = asin(clamp(z / r, -1.0, 1.0))
    lon = atan(y, x)

    return lat, lon
end
