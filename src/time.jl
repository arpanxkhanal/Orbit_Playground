function mean_motion(a::Real, μ::Real = EARTH_MU)
    a > 0 || throw(ArgumentError("Semi-major axis must be positive."))
    μ > 0 || throw(ArgumentError("Gravitational parameter must be positive."))
    return sqrt(μ / a^3)
end

function orbital_period(a::Real, μ::Real = EARTH_MU)
    return 2π / mean_motion(a, μ)
end

function mean_anomaly(t::Real, n::Real, M0::Real = 0.0)
    return M0 + n * t
end

function gast(t::Real, ωE::Real = EARTH_ROTATION_RATE)
    return ωE * t
end
