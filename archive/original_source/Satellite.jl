"""
Satellite Structure

Stores the Keplerian elements of a satellite.

Fields:
name : Satellite name
a    : Semi-major axis [km]
e    : Eccentricity [-]
i    : Inclination [rad]
Ω    : Right Ascension of the Ascending Node [rad]
ω    : Argument of Perigee [rad]
M0   : Mean anomaly at epoch [rad]
"""

struct Satellite

    name::String

    a::Float64

    e::Float64

    i::Float64

    Ω::Float64

    ω::Float64

    M0::Float64

end

struct SatelliteState

    satellite::Satellite

    position_eci::Vector{Float64}

    velocity_eci::Vector{Float64}

    position_ecef::Vector{Float64}

    latitude::Float64

    longitude::Float64
    
end