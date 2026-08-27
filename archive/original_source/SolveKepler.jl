# M = E - esinE
# Delta M = Delta E (1 - e Cos E)
# E (first) = M + e sin(M) + e2/2 sin(2M)

function solve_kepler(e, M)
    
    if e < 0 || e >= 1
        error("Eccentricity must be in the range [0, 1)")

    end


    E = M + e * sin(M) + e^2/2 * sin(2M)
    tol = 1e-12

    for k in 1:50 

        dM = M - (E - e * sin(E) ) 
        dE = dM / (1 - e * cos(E))
        E = E + dE
        if abs(dE) < tol
            return E
        end
    end
    error("Kepler's equation did not converge after 50 iterations")
end