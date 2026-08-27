function closure_error_matrix(
    sat::Satellite,
    repeat_period_value::Real;
    μ::Real = EARTH_MU,
    ωE::Real = EARTH_ROTATION_RATE,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
    n_cycles::Int = 5,
    n_samples::Int = 200,
)
    repeat_period_value > 0 || throw(ArgumentError("repeat_period must be positive."))
    n_cycles > 0 || throw(ArgumentError("n_cycles must be positive."))
    n_samples >= 2 || throw(ArgumentError("n_samples must be at least 2."))

    t_reference = collect(range(
        0.0,
        repeat_period_value;
        length = n_samples,
    ))

    reference = Matrix{Float64}(undef, 3, n_samples)
    errors = Matrix{Float64}(undef, n_cycles, n_samples)

    for (j, t) in enumerate(t_reference)
        reference[:, j] = compute_state(
            sat, t, μ, ωE;
            use_j2 = true,
            Re = Re,
            J2 = J2,
        ).position_ecef
    end

    for cycle in 1:n_cycles
        for (j, t) in enumerate(t_reference)
            current = compute_state(
                sat,
                t + cycle * repeat_period_value,
                μ,
                ωE;
                use_j2 = true,
                Re = Re,
                J2 = J2,
            ).position_ecef

            errors[cycle, j] = norm(current - reference[:, j])
        end
    end

    return errors
end

function closure_statistics(errors::AbstractMatrix)
    rms = sqrt.(vec(mean(errors .^ 2, dims = 2)))
    maximum_error = vec(maximum(errors, dims = 2))
    mean_error = vec(mean(errors, dims = 2))

    return (
        mean = mean_error,
        rms = rms,
        maximum = maximum_error,
    )
end

function repeat_cycle_closure(
    sat::Satellite,
    repeat_period_value::Real;
    kwargs...,
)
    errors = closure_error_matrix(
        sat,
        repeat_period_value;
        kwargs...,
    )

    stats = closure_statistics(errors)

    return (
        errors = errors,
        mean_error = stats.mean,
        rms_error = stats.rms,
        maximum_error = stats.maximum,
    )
end

function element_rate_check(
    sat::Satellite,
    t::Real;
    μ::Real = EARTH_MU,
    Re::Real = EARTH_RADIUS,
    J2::Real = EARTH_J2,
)
    Ωdot, ωdot, Mdot = compute_j2_secular_rates(
        sat.a, sat.e, sat.i;
        μ = μ, Re = Re, J2 = J2,
    )

    propagated = propagate_j2_elements(
        sat, t;
        μ = μ, Re = Re, J2 = J2,
    )

    predicted = [
        sat.Ω + Ωdot * t,
        sat.ω + ωdot * t,
        sat.M0 + Mdot * t,
    ]

    actual = [propagated.Ω, propagated.ω, propagated.M0]

    return (
        predicted = predicted,
        actual = actual,
        residual = actual .- predicted,
        rates = (Ωdot = Ωdot, ωdot = ωdot, Mdot = Mdot),
    )
end
