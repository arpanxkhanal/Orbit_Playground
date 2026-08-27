module SatelliteOrbitPropagation

using LinearAlgebra
using Statistics
using Plots

include("constants.jl")
include("types.jl")
include("time.jl")
include("rotations.jl")
include("kepler.jl")
include("coordinates.jl")
include("perturbations.jl")
include("propagation.jl")
include("repeat_orbit.jl")
include("analysis.jl")
include("plotting.jl")
include("validation.jl")

export
    # Constants
    EARTH_MU, EARTH_RADIUS, EARTH_J2, EARTH_ROTATION_RATE,

    # Types
    Satellite, SatelliteState, OrbitalElementHistory, PropagationResult,

    # Time / Kepler
    mean_motion, orbital_period, mean_anomaly, gast,
    solve_kepler,

    # Coordinate transformations
    kepler_to_pqw, pqw_to_eci, eci_to_ecef, ecef_to_latlon, R3,

    # Perturbations
    compute_j2_secular_rates, propagate_j2_elements,

    # State / propagation
    compute_state, propagate_kepler, propagate_j2,

    # Repeat-orbit design
    repeat_ratio, repeat_orbit_error, design_repeat_orbit,
    repeat_period,

    # Analysis
    closure_error_matrix, closure_statistics, repeat_cycle_closure,
    element_rate_check,

    # Plotting
    plot_orbit, plot_groundtrack, plot_orbital_elements,

    # Reporting / validation
    summary, validation_status

end
