"""
    summary(result)

Prints a summary of a propagated orbit.

Inputs:
    result : PropagationResult

Returns:
    Nothing
"""
function summary(result::PropagationResult)

    println()
    println("================== PROPAGATION SUMMARY ==================")
    println()

    println("Satellite              : ", result.satellite)
    println("Propagator             : ", result.propagator)
    println("Time Step              : ", result.dt, " s")
    println("Revolutions            : ", result.revolutions)
    println(
        "Propagation Time       : ",
        round(result.time[end], digits = 2),
        " s (",
        round(result.time[end] / 3600, digits = 2),
        " hr)"
    )
    println("Stored Epochs          : ", length(result.time))

    println()

    println("---------------------------------------------------------")
    println("Initial Orbital Elements")
    println("---------------------------------------------------------")
    println("a = $(result.elements.a[1]) km")
    println("e = $(result.elements.e[1])")
    println("i = $(round(rad2deg(result.elements.i[1]), digits=6))°")
    println("Ω = $(round(rad2deg(result.elements.Ω[1]), digits=6))°")
    println("ω = $(round(rad2deg(result.elements.ω[1]), digits=6))°")
    println("M = $(round(rad2deg(result.elements.M[1]), digits=6))°")

    println()

    println("---------------------------------------------------------")    
    println("Final Orbital Elements")
    println("---------------------------------------------------------")
    println("a = $(result.elements.a[end]) km")
    println("e = $(result.elements.e[end])")
    println("i = $(round(rad2deg(result.elements.i[end]), digits=6))°")
    println("Ω = $(round(rad2deg(result.elements.Ω[end]), digits=6))°")
    println("ω = $(round(rad2deg(result.elements.ω[end]), digits=6))°")
    println("M = $(round(rad2deg(result.elements.M[end]), digits=6))°")

    println()

    println("---------------------------------------------------------")
    println("Element Changes")
    println("---------------------------------------------------------")
    println("Δa = $(result.elements.a[end] - result.elements.a[1]) km")
    println("Δe = $(result.elements.e[end] - result.elements.e[1])")
    println("Δi = $(round(rad2deg(result.elements.i[end] - result.elements.i[1]), digits=6))°")
    println("ΔΩ = $(round(rad2deg(result.elements.Ω[end] - result.elements.Ω[1]), digits=6))°")
    println("Δω = $(round(rad2deg(result.elements.ω[end] - result.elements.ω[1]), digits=6))°")
    println("ΔM = $(round(rad2deg(result.elements.M[end] - result.elements.M[1]), digits=6))°")

    println("=========================================================")

    return nothing

end
