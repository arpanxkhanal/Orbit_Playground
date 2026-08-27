function _time_grid(T::Real, dt::Real, revolutions::Int)
    dt > 0 || throw(ArgumentError("dt must be positive."))
    revolutions > 0 || throw(ArgumentError("revolutions must be positive."))

    t_end = revolutions * T
    t = collect(0.0:dt:t_end)

    # Ensure the requested endpoint is represented exactly.
    if isempty(t) || !isapprox(t[end], t_end; atol = 10eps(Float64) * max(1.0, t_end))
        push!(t, t_end)
    else
        t[end] = t_end
    end

    return t
end

function _build_result(
    propagator::String,
    sat::Satellite,
    t::Vector{Float64},
    states::Vector{SatelliteState},
    dt::Float64,
    revolutions::Int,
    repeat_revolutions::Int,
)
    N = length(t)

    positions_eci = Matrix{Float64}(undef, 3, N)
    velocities_eci = Matrix{Float64}(undef, 3, N)
    positions_ecef = Matrix{Float64}(undef, 3, N)
    latitudes = Vector{Float64}(undef, N)
    longitudes = Vector{Float64}(undef, N)

    a_history = Vector{Float64}(undef, N)
    e_history = Vector{Float64}(undef, N)
    i_history = Vector{Float64}(undef, N)
    Ω_history = Vector{Float64}(undef, N)
    ω_history = Vector{Float64}(undef, N)
    M_history = Vector{Float64}(undef, N)

    for k in eachindex(states)
        s = states[k]
        positions_eci[:, k] = s.position_eci
        velocities_eci[:, k] = s.velocity_eci
        positions_ecef[:, k] = s.position_ecef
        latitudes[k] = s.latitude
        longitudes[k] = s.longitude

        a_history[k] = s.satellite.a
        e_history[k] = s.satellite.e
        i_history[k] = s.satellite.i
        Ω_history[k] = s.satellite.Ω
        ω_history[k] = s.satellite.ω
        M_history[k] = s.satellite.M0
    end

    elements = OrbitalElementHistory(
        a_history,
        e_history,
        i_history,
        Ω_history,
        ω_history,
        M_history,
    )

    return PropagationResult(
        propagator,
        dt,
        revolutions,
        sat.name,
        repeat_revolutions,
        t,
        positions_eci,
        velocities_eci,
        positions_ecef,
        latitudes,
        longitudes,
        elements,
    )
end

function compute_state(
    sat::Satellite,
    t::Real,
    μ::Real = EARTH_MU,
    ωE::Real = EARTH_ROTATION_RATE;
    use_j2::Bool = true,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
)
    sat_now = use_j2 ?
        propagate_j2_elements(sat, t; μ = μ, Re = Re, J2 = J2) :
        sat

    n = mean_motion(sat.a, μ)
    M = use_j2 ? sat_now.M0 : mean_anomaly(t, n, sat.M0)
    E = solve_kepler(sat_now.e, M)

    r_pqw, v_pqw = kepler_to_pqw(sat_now, E, μ)
    r_eci, v_eci = pqw_to_eci(sat_now, r_pqw, v_pqw)

    r_ecef = eci_to_ecef(r_eci, gast(t, ωE))
    lat, lon = ecef_to_latlon(r_ecef)

    return SatelliteState(
        sat_now,
        r_eci,
        v_eci,
        r_ecef,
        rad2deg(lat),
        rad2deg(lon),
    )
end

function propagate_kepler(
    sat::Satellite;
    μ::Real = EARTH_MU,
    ωE::Real = EARTH_ROTATION_RATE,
    dt::Real = 10.0,
    revolutions::Int = 2,
)
    T = orbital_period(sat.a, μ)
    t = _time_grid(T, dt, revolutions)

    states = [
        compute_state(sat, ti, μ, ωE; use_j2 = false)
        for ti in t
    ]

    return _build_result(
        "Analytical Keplerian",
        sat,
        t,
        states,
        Float64(dt),
        revolutions,
        revolutions,
    )
end

function propagate_j2(
    sat::Satellite;
    μ::Real = EARTH_MU,
    ωE::Real = EARTH_ROTATION_RATE,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
    dt::Real = 10.0,
    revolutions::Int = 2,
    repeat_revolutions::Int = revolutions,
)
    T = orbital_period(sat.a, μ)
    t = _time_grid(T, dt, revolutions)

    states = [
        compute_state(
            sat, ti, μ, ωE;
            use_j2 = true,
            Re = Re,
            J2 = J2,
        )
        for ti in t
    ]

    return _build_result(
        "Analytical secular J2",
        sat,
        t,
        states,
        Float64(dt),
        revolutions,
        repeat_revolutions,
    )
end
