Base.@kwdef struct Satellite
    name::String
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    M0::Float64
end

struct SatelliteState
    satellite::Satellite
    position_eci::Vector{Float64}
    velocity_eci::Vector{Float64}
    position_ecef::Vector{Float64}
    latitude::Float64
    longitude::Float64
end

struct OrbitalElementHistory
    a::Vector{Float64}
    e::Vector{Float64}
    i::Vector{Float64}
    Ω::Vector{Float64}
    ω::Vector{Float64}
    M::Vector{Float64}
end

struct PropagationResult
    propagator::String
    dt::Float64
    revolutions::Int
    satellite::String
    repeat_revolutions::Int
    time::Vector{Float64}
    positions_eci::Matrix{Float64}
    velocities_eci::Matrix{Float64}
    positions_ecef::Matrix{Float64}
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}
    elements::OrbitalElementHistory
end
