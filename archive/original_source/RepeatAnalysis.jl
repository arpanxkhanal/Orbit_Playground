using LinearAlgebra
using Statistics

function closure_error(result::PropagationResult)

    cycle_number = Int[]
    rms_errors = Float64[]
    max_errors = Float64[]
    
    # Orbital period
    T = result.time[end] / result.revolutions

    # Samples
    samples_per_rev = round(Int, T / result.dt)
    samples_per_cycle = result.repeat_revolutions * samples_per_rev + 1

    # Number of complete repeat cycles
    num_cycles = div(result.revolutions, result.repeat_revolutions)

    # Reference trajectory (always the first repeat cycle)
    reference = result.positions_ecef[:, 1:samples_per_cycle]

    println("Orbital Period = ", T)
    println("Samples per Revolution = ", samples_per_rev)
    println("Samples per Cycle = ", samples_per_cycle)
    println("Number of Repeat Cycles = ", num_cycles)
    println("Reference Size = ", size(reference))

    # Loop over every repeat cycle after the reference
    for cycle = 2:num_cycles

        start_idx = 1 + (cycle - 1) * (samples_per_cycle - 1)
        stop_idx  = start_idx + samples_per_cycle - 1

        current = result.positions_ecef[:, start_idx:stop_idx]

        distances = zeros(Float64, samples_per_cycle)

        for i = 1:samples_per_cycle
            distances[i] = norm(reference[:, i] - current[:, i])
        end

        rms_error = sqrt(mean(distances .^ 2))
        max_error = maximum(distances)
        push!(cycle_number, cycle)
        push!(rms_errors, rms_error)
        push!(max_errors, max_error)

        println("\n------------------------------")
        println("Cycle ", cycle)
        println("------------------------------")
        println("Start Index = ", start_idx)
        println("Stop Index  = ", stop_idx)
        println("Current Size = ", size(current))

        println("Start Time = ", result.time[start_idx])
        println("Stop Time  = ", result.time[stop_idx])
        println("Distances = ", distances)
        println("RMS Closure Error = ", rms_error)
        println("Maximum Closure Error = ", max_error)

    end
    return cycle_number, rms_errors, max_errors

end