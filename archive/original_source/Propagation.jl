"""
function propagate_kepler(sat, GM, ωE, dt, revolutions)

Propagates a satellite using the analytical two-body (Keplerian)
orbit model.

Inputs:
sat         : Satellite structure
GM          : Gravitational parameter [km^3/s^2]
ωE          : Earth's rotation rate [rad/s]
dt          : Time step [s]
revolutions : Number of orbital revolutions

Returns:
positions_eci  : 3×N matrix of ECI position vectors [km]
positions_ecef : 3×N matrix of ECEF position vectors [km]
velocities_eci : 3×N matrix of ECI velocity vectors [km/s]
latitudes      : Latitude vector [deg]
longitudes     : Longitude vector [deg]

"""
function propagate_kepler(
    sat::Satellite,
    GM::Float64,
    ωE::Float64;
    dt::Float64 = 10.0,
    revolutions::Int = 2
)
    # Compute Mean Motion and ORbital Peroid
    n = mean_motion(sat.a, GM)
    T = orbital_period(sat.a, GM)

    # Time vector
    t_list = collect(0.0:dt:revolutions*T)

    # Number of time steps
    N = length(t_list)

    # Preallocate memory
    positions_eci = Matrix{Float64}(undef, 3, N)
    positions_ecef = Matrix{Float64}(undef, 3, N)
    velocities_eci = Matrix{Float64}(undef, 3, N)
    latitudes = Vector{Float64}(undef, N)
    longitudes = Vector{Float64}(undef, N) 

    # Propagation loop
    for (k, t) in enumerate(t_list)

        M = mean_anomaly(t, n, sat.M0)
        E = solve_kepler(sat.e, M)

        r_pqw, v_pqw = kepler_to_pqw(sat, E, GM)

        r_eci, v_eci = pqw_to_eci(sat, r_pqw, v_pqw)

        GAST = gast(t, ωE)

        r_ecef = eci_to_ecef(r_eci, GAST)

        lat, lon = ecef_to_latlon(r_ecef)

        # Store Results
        positions_eci[:, k] = r_eci
        velocities_eci[:, k] = v_eci
        positions_ecef[:, k] = r_ecef
        latitudes[k] = rad2deg(lat)
        longitudes[k] = rad2deg(lon)

    end

    return PropagationResult(
        t_list,
        positions_eci,
        velocities_eci,
        positions_ecef,
        latitudes,
        longitudes
    )




end



