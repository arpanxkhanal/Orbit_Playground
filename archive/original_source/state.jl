"""
State.jl

This function is aimed to give me the 
eci ecef position and lat long for any given epoch t

compute_state(sat, t, GM, ωE)

keplerian elements will be propageted using the J2 perturbation
model. Then conversion to eci, ecef and lat long will be done.

Inputs:
sat : Satellite structure
t   : Time since epoch [s]
GM  : Gravitational parameter [km³/s²]
ωE  : Earth's rotation rate [rad/s]

Returns:
r_eci : ECI position vector [km]
v_eci : ECI velocity vector [km/s]
r_ecef : ECEF position vector [km]
lat : Latitude [deg]
lon : Longitude [deg]

"""

function compute_state(
    sat::Satellite,
    t::Float64,
    GM::Float64,
    ωE::Float64 
)

    # Propagate the orbital elements using J2 secular model
    sat_now = propagate_j2_elements(sat, t)

    # Solve Kepler's equation
    E = solve_kepler(sat_now.e, sat_now.M0)

    # Compute the satellite state in PQW frame 
    r_pqw, v_pqw = kepler_to_pqw(sat_now, E, GM)

    # Convert to ECI frame
    r_eci, v_eci = pqw_to_eci(sat_now, r_pqw, v_pqw)

    # Convert to ECEF frame
    GAST = gast(t, ωE)
    r_ecef = eci_to_ecef(r_eci, GAST)

    # Convert to Latitude and Longitude
    lat, lon = ecef_to_latlon(r_ecef) 

    return SatelliteState(
        sat_now,
        r_eci,
        v_eci,
        r_ecef,
        rad2deg(lat),
        rad2deg(lon)
    )
end