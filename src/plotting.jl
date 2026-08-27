function _break_groundtrack(lat, lon)
    lat_plot = Float64[lat[1]]
    lon_plot = Float64[lon[1]]

    for k in 2:length(lat)
        Δlon = abs(lon[k] - lon[k - 1])
        wrap = Δlon > 180
        pole = (abs(lat[k]) > 88 || abs(lat[k - 1]) > 88) && Δlon > 20

        if wrap || pole
            push!(lat_plot, NaN)
            push!(lon_plot, NaN)
        end

        push!(lat_plot, lat[k])
        push!(lon_plot, lon[k])
    end

    return lat_plot, lon_plot
end

function plot_orbit(
    positions::AbstractMatrix;
    title::String = "Orbit",
)
    size(positions, 1) == 3 || throw(ArgumentError(
        "positions must be a 3×N matrix."
    ))

    p = plot3d(
        positions[1, :],
        positions[2, :],
        positions[3, :],
        xlabel = "X (km)",
        ylabel = "Y (km)",
        zlabel = "Z (km)",
        title = title,
        linewidth = 2,
        legend = false,
    )

    scatter!(
        p,
        [positions[1, 1]],
        [positions[2, 1]],
        [positions[3, 1]],
        label = "Start",
        markersize = 4,
    )

    return p
end

function plot_groundtrack(
    latitudes::AbstractVector,
    longitudes::AbstractVector;
    title::String = "Ground Track",
)
    length(latitudes) == length(longitudes) || throw(ArgumentError(
        "Latitude and longitude arrays must have equal length."
    ))

    lat_plot, lon_plot = _break_groundtrack(latitudes, longitudes)

    return plot(
        lon_plot,
        lat_plot,
        xlabel = "Longitude (°)",
        ylabel = "Latitude (°)",
        title = title,
        xlims = (-180, 180),
        ylims = (-90, 90),
        xticks = -180:60:180,
        yticks = -90:30:90,
        aspect_ratio = :equal,
        linewidth = 1.5,
        legend = false,
        grid = false,
        framestyle = :box,
    )
end

function plot_groundtrack(result::PropagationResult; title::String = "Ground Track")
    return plot_groundtrack(
        result.latitudes,
        result.longitudes;
        title = title,
    )
end

function plot_orbital_elements(result::PropagationResult)
    t = result.time ./ 3600
    e = result.elements

    return [
        plot(t, e.a, xlabel = "Time (hr)", ylabel = "a (km)",
             title = "Semi-major Axis", legend = false),
        plot(t, e.e, xlabel = "Time (hr)", ylabel = "e",
             title = "Eccentricity", legend = false),
        plot(t, rad2deg.(e.i), xlabel = "Time (hr)", ylabel = "i (deg)",
             title = "Inclination", legend = false),
        plot(t, rad2deg.(e.Ω), xlabel = "Time (hr)", ylabel = "Ω (deg)",
             title = "RAAN", legend = false),
        plot(t, rad2deg.(e.ω), xlabel = "Time (hr)", ylabel = "ω (deg)",
             title = "Argument of Perigee", legend = false),
        plot(t, rad2deg.(e.M), xlabel = "Time (hr)", ylabel = "M (deg)",
             title = "Mean Anomaly", legend = false),
    ]
end

function summary(result::PropagationResult)
    e = result.elements

    println("========== PROPAGATION SUMMARY ==========")
    println("Satellite       : ", result.satellite)
    println("Propagator      : ", result.propagator)
    println("Time step       : ", result.dt, " s")
    println("Revolutions     : ", result.revolutions)
    println("Duration        : ", round(result.time[end], digits = 3), " s")
    println("Stored epochs   : ", length(result.time))
    println()
    println("Initial elements")
    println("a = ", e.a[1], " km")
    println("e = ", e.e[1])
    println("i = ", rad2deg(e.i[1]), " deg")
    println("Ω = ", rad2deg(e.Ω[1]), " deg")
    println("ω = ", rad2deg(e.ω[1]), " deg")
    println("M = ", rad2deg(e.M[1]), " deg")
    println()
    println("Final elements")
    println("a = ", e.a[end], " km")
    println("e = ", e.e[end])
    println("i = ", rad2deg(e.i[end]), " deg")
    println("Ω = ", rad2deg(e.Ω[end]), " deg")
    println("ω = ", rad2deg(e.ω[end]), " deg")
    println("M = ", rad2deg(e.M[end]), " deg")
    println("==========================================")

    return nothing
end
