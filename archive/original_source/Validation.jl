"""
compute_closure_errors(
    sat,
    GM,
    ωE,
    repeat_period;
    n_cycles = 10,
    n_samples = 100
)

Computes the repeat-orbit closure error over multiple repeat cycles.

The repeat period is divided into equally spaced sample epochs.
For each sample epoch, the reference position is computed once.
The satellite position is then recomputed after successive repeat
periods and compared against the reference position.

Inputs:
sat           : Satellite structure
GM            : Earth's gravitational parameter [km³/s²]
ωE            : Earth's rotation rate [rad/s]
repeat_period : Repeat period [s]

Keyword Arguments:
n_cycles  : Number of repeat cycles to evaluate
n_samples : Number of sample points within one repeat period

Returns:
errors : Matrix of closure errors [km]

Rows    -> Repeat cycle
Columns -> Sample location within repeat period
"""
function compute_closure_errors(
    sat::Satellite,
    GM::Float64,
    ωE::Float64,
    repeat_period::Float64;
    n_cycles::Int = 10,
    n_samples::Int = 100
)

    # --------------------------------------------------
    # Reference sample times
    # --------------------------------------------------

    t_reference = collect(
        range(
            0.0,
            repeat_period;
            length = n_samples
        )
    )

    # --------------------------------------------------
    # Allocate memory
    # --------------------------------------------------

    reference_positions = Matrix{Float64}(undef, 3, n_samples)

    errors = Matrix{Float64}(undef, n_cycles, n_samples)

    # --------------------------------------------------
    # Compute reference positions
    # --------------------------------------------------

    for (i,t) in enumerate(t_reference)

        state = compute_state(
            sat,
            t,
            GM,
            ωE
        )

        reference_positions[:,i] = state.position_ecef

    end

    # --------------------------------------------------
    # Compute closure errors
    # --------------------------------------------------

    for cycle in 1:n_cycles

        for (i,t) in enumerate(t_reference)

            comparison_time = t + cycle * repeat_period

            state = compute_state(
                sat,
                comparison_time,
                GM,
                ωE
            )

            errors[cycle,i] = norm(
                state.position_ecef -
                reference_positions[:,i]
            )

        end

    end

    return errors

end

function closure_statistics(errors)

    mean_error = vec(mean(errors, dims=2))

    std_error = vec(std(errors, dims=2))

    return mean_error, std_error

end