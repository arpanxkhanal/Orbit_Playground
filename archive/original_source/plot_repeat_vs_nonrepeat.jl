
"""
Breaks the plotted line whenever the ground track crosses
the International Date Line or passes over the poles.
"""
function break_groundtrack(lat, lon)

    lat_new = Float64[]
    lon_new = Float64[]

    push!(lat_new, lat[1])
    push!(lon_new, lon[1])

    for k in 2:length(lat)

        Δlon = abs(lon[k] - lon[k-1])

        # Break at date line
        wrap = Δlon > 180

        # Break near poles
        pole = (abs(lat[k]) > 88 || abs(lat[k-1]) > 88) && Δlon > 20

        if wrap || pole
            push!(lat_new, NaN)
            push!(lon_new, NaN)
        end

        push!(lat_new, lat[k])
        push!(lon_new, lon[k])

    end

    return lat_new, lon_new

end


function plot_revolutions(
    result::PropagationResult,
    T;
    title = "Ground Tracks"
)

    lat_plot, lon_plot = break_groundtrack(
        result.latitudes,
        result.longitudes
    )

    p = plot(
        lon_plot,
        lat_plot,

        title = title,

        xlabel = "Longitude (°)",
        ylabel = "Latitude (°)",

        xlims = (-180,180),
        ylims = (-90,90),

        xticks = -180:60:180,
        yticks = -90:30:90,

        aspect_ratio = :equal,

        legend = false,
        grid = false,

        framestyle = :box,

        background_color = :white,

        linewidth = 1.2,

        color = :royalblue
    )

    return p

end