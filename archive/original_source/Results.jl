

"""
OrbitalElementHistory

Stores the orbital elements throughout propagation.

Fields:
a : Semi-major axis [km]
e : Eccentricity [-]
i : Inclination [rad]
Ω : Right Ascension of Ascending Node [rad]
ω : Argument of Perigee [rad]
M : Mean Anomaly [rad]
"""
struct OrbitalElementHistory

    a::Vector{Float64}

    e::Vector{Float64}

    i::Vector{Float64}

    Ω::Vector{Float64}

    ω::Vector{Float64}

    M::Vector{Float64}

end



"""

PropagationResult

Stores the propagated state vectors and derived quantities.

Fields:
time            : Time vector [s]
positions_eci   : 3×N ECI position matrix [km]
velocities_eci  : 3×N ECI velocity matrix [km/s]
positions_ecef  : 3×N ECEF position matrix [km]
latitudes       : Latitude vector [deg]
longitudes      : Longitude vector [deg]
elements        : OrbitalElementHistory
"""

struct PropagationResult

    # Metadata
    propagator::String
    dt::Float64
    revolutions::Int
    satellite::String
    repeat_revolutions::Int

    # Time
    time::Vector{Float64}

    # States
    positions_eci::Matrix{Float64}
    velocities_eci::Matrix{Float64}
    positions_ecef::Matrix{Float64}

    # Ground Track
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}

    # Orbital Elements
    elements::OrbitalElementHistory

end