"""
function propagate_j2(sat, GM, ωE, dt, revolutions)

Propagates a satellite using the analytical secular J₂
perturbation model.

Inputs:
sat         : Satellite structure
GM          : Gravitational parameter [km³/s²]
ωE          : Earth's rotation rate [rad/s]
dt          : Time step [s]
revolutions : Number of orbital revolutions

Returns:
PropagationResult containing

- Time history
- ECI and ECEF state vectors
- Ground track (latitude and longitude)
- Orbital element history

"""
function propagate_j2_old(
    sat::Satellite,
    GM::Float64,
    ωE::Float64;
    dt::Float64 = 10.0,
    revolutions::Int = 2,
    repeat_revolutions::Int = revolutions
)
    # Orbital Period
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
    a_history = Vector{Float64}(undef, N)
    e_history = Vector{Float64}(undef, N)
    i_history = Vector{Float64}(undef, N)
    Ω_history = Vector{Float64}(undef, N)
    ω_history = Vector{Float64}(undef, N)
    M_history = Vector{Float64}(undef, N)

    # Propagation loop
    for (k,t) in enumerate(t_list)

        # Propagate the orbital elements using the J2 secular model
        sat_now = propagate_j2_elements(sat,t)



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

        # Store Results
        positions_eci[:, k] = r_eci
        velocities_eci[:, k] = v_eci
        positions_ecef[:, k] = r_ecef
        latitudes[k] = rad2deg(lat)
        longitudes[k] = rad2deg(lon)
        
        # Orbital Element History

        a_history[k] = sat_now.a
        e_history[k] = sat_now.e
        i_history[k] = sat_now.i
        Ω_history[k] = sat_now.Ω
        ω_history[k] = sat_now.ω
        M_history[k] = sat_now.M0

    end

    elements = OrbitalElementHistory(
        a_history,
        e_history,
        i_history,
        Ω_history,
        ω_history,
        M_history
    )


    return PropagationResult(
    "Analytical J2",
    dt,
    revolutions,
    sat.name,
    repeat_revolutions,
    t_list,
    positions_eci,
    velocities_eci,
    positions_ecef,
    latitudes,
    longitudes,
    elements
    )
end
