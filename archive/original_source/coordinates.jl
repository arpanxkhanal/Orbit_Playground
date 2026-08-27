""" 
function kepler_to_pqw(sat, E, GM) 
This takes in : 

    1. sat : Satellite structure which contains the orbital elements,
        name : Satellite name
        a    : Semi-major axis [km]
        e    : Eccentricity [-]
        i    : Inclination [rad]
        Ω    : Right Ascension of the Ascending Node [rad]
        ω    : Argument of Perigee [rad]
        M0   : Mean anomaly at epoch [rad]
    
    2. E : Eccentric anomaly [rad]

    3. GM : Gravitational parameter [km^3/s^2]

Returns :

    1. r_pqw : position vector in PQW frame [km]
    2. v_pqw : velocity vector in PQW frame [km/s]

end




function pqw_to_eci(sat, r_pqw, v_pqw)
    This takes in :
    1. sat : Satellite structure which contains the orbital elements,
        name : Satellite name
        a    : Semi-major axis [km]
        e    : Eccentricity [-]
        i    : Inclination [rad]
        Ω    : Right Ascension of the Ascending Node [rad]
        ω    : Argument of Perigee [rad]
        M0   : Mean anomaly at epoch [rad]

    2. r_pqw : position vector in PQW frame [km]

    3. v_pqw : velocity vector in PQW frame [km/s]

    Returns :

    1. r_eci : position vector in ECI frame [km]

    2. v_eci : velocity vector in ECI frame [km/s]

end


function eci_to_ecef(r_eci, GAST)

Converts the position vector from the Earth-Centered Inertial (ECI)
frame to the Earth-Centered Earth-Fixed (ECEF) frame.

Inputs:
r_eci : Position vector in ECI [km]
GAST  : Greenwich Apparent Sidereal Time [rad]

Returns:
r_ecef : Position vector in ECEF [km]


function ecef_to_latlon(r_ecef)

Converts an Earth-Centered Earth-Fixed (ECEF) position vector
to geocentric latitude and longitude assuming a spherical Earth.

Inputs:
r_ecef : Position vector in ECEF [km]

Returns:
lat : Latitude [rad]
lon : Longitude [rad]



"""

using LinearAlgebra

function kepler_to_pqw(sat, E, GM)

    a = sat.a
    e = sat.e

    x = a * (cos(E) - e)
    y = a * sqrt(1 - e^2) * sin(E)
    z = 0.0

    r_pqw = [x, y, z]

    r = norm(r_pqw)

    f = sqrt(GM * a) / r

    vx = -f * sin(E)
    vy =  f * sqrt(1 - e^2) * cos(E)
    vz = 0.0

    v_pqw = [vx, vy, vz]

    return r_pqw, v_pqw

end



function pqw_to_eci(sat, r_pqw, v_pqw)

    i = sat.i
    Ω = sat.Ω
    ω = sat.ω

    Q = R3(-Ω) * R1(-i) * R3(-ω)

    r_eci = Q * r_pqw

    v_eci = Q * v_pqw

    return r_eci, v_eci

end



function eci_to_ecef(r_eci, GAST)

    r_ecef = R3(GAST) * r_eci

    return r_ecef

end


function ecef_to_latlon(r_ecef)

    x, y, z = r_ecef

    r = norm(r_ecef)

    lat = asin(z / r)

    lon = atan(y, x)

    return lat, lon

end