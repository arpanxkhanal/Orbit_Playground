function R1(theta)
    return [
        1.0 0.0 0.0 
        0.0 cos(theta) sin(theta)
        0.0 -sin(theta) cos(theta)
    ]
end

function R3(theta)
    return [
        cos(theta) sin(theta) 0.0
        -sin(theta) cos(theta) 0.0
        0.0 0.0 1.0
    ]
end
