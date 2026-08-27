using Plots

"""
plot_orbit(positions; title="Orbit")

Plots a 3D satellite orbit.

Inputs:
positions : 3×N position matrix [km]

Keyword Arguments:
title : Plot title

Returns:
3D orbit plot
"""
function plot_orbit(
    positions;
    title="Orbit"
)

    x = positions[1,:]
    y = positions[2,:]
    z = positions[3,:]

    p = plot3d(
        x,
        y,
        z,
        xlabel = "X (km)",
        ylabel = "Y (km)",
        zlabel = "Z (km)",
        linewidth = 2,
        legend = false,
        title = title
    )

    scatter!(
        p,
        [x[1]],
        [y[1]],
        [z[1]],
        markersize = 4,
        label = "Start"
    )

    display(p)

    return p

end

"""
plot_groundtrack(latitudes, longitudes)

Plots the satellite ground track.

Inputs:
latitudes  : Latitude vector [deg]
longitudes : Longitude vector [deg]

Returns:
Ground track plot
"""
function plot_groundtrack(
    latitudes,
    longitudes;
    title="Ground Track"
)
    # Break the line when crossing the ±180° meridian

    plot_lon = Float64[]
    plot_lat = Float64[]

    push!(plot_lon, longitudes[1])
    push!(plot_lat, latitudes[1])

    for k in 2:length(longitudes)

        if abs(longitudes[k] - longitudes[k-1]) > 180
            push!(plot_lon, NaN)
            push!(plot_lat, NaN)
        end

        push!(plot_lon, longitudes[k])
        push!(plot_lat, latitudes[k])

    end
    p = plot(
        plot_lon,
        plot_lat,
        xlabel = "Longitude (deg)",
        ylabel = "Latitude (deg)",
        linewidth = 2,
        legend = false,
        title = title
    )

    scatter!(
        p,
        [longitudes[1]],
        [latitudes[1]],
        markersize = 4,
        label = "Start"
    )

    xlims!(-180,180)

    ylims!(-90,90)

    display(p)

    return p

end