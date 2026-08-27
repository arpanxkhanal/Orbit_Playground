function R1(θ::Real)
    return [
        1.0 0.0 0.0
        0.0 cos(θ) sin(θ)
        0.0 -sin(θ) cos(θ)
    ]
end

function R3(θ::Real)
    return [
        cos(θ) sin(θ) 0.0
        -sin(θ) cos(θ) 0.0
        0.0 0.0 1.0
    ]
end
