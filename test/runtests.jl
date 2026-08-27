using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "SatelliteOrbitPropagation.jl"))
using .SatelliteOrbitPropagation

@testset "Satellite Orbit Propagation" begin
    sat = Satellite(
        "TestSat",
        7000.0,
        0.01,
        deg2rad(98.0),
        deg2rad(10.0),
        deg2rad(20.0),
        deg2rad(30.0),
    )

    @testset "Kepler solver" begin
        e = 0.1
        M = 1.2
        E = solve_kepler(e, M)
        @test isapprox(E - e * sin(E), M; atol = 1e-11)
        @test_throws ArgumentError solve_kepler(-0.1, M)
        @test_throws ArgumentError solve_kepler(1.0, M)
    end

    @testset "Time relations" begin
        n = mean_motion(sat.a)
        T = orbital_period(sat.a)
        @test isapprox(n * T, 2π; atol = 1e-12)
        @test mean_anomaly(10.0, n, sat.M0) == sat.M0 + 10n
    end

    @testset "Rotations" begin
        Q = R3(0.37)
        @test isapprox(Q' * Q, Matrix(I, 3, 3); atol = 1e-12)
        @test isapprox(det(Q), 1.0; atol = 1e-12)
    end

    @testset "Coordinate/state construction" begin
        r_pqw, v_pqw = kepler_to_pqw(sat, 0.4)
        @test length(r_pqw) == 3
        @test length(v_pqw) == 3

        state = compute_state(sat, 0.0; use_j2 = true)
        @test size(state.position_eci) == (3,)
        @test size(state.position_ecef) == (3,)
        @test -90 <= state.latitude <= 90
        @test -180 <= state.longitude <= 180
    end

    @testset "Propagation result" begin
        result = propagate_kepler(sat; dt = 60.0, revolutions = 2)
        @test length(result.time) > 2
        @test size(result.positions_eci, 1) == 3
        @test size(result.positions_ecef, 1) == 3
        @test length(result.latitudes) == length(result.time)
        @test length(result.elements.Ω) == length(result.time)

        result_j2 = propagate_j2(
            sat;
            dt = 60.0,
            revolutions = 2,
            repeat_revolutions = 2,
        )
        @test result_j2.propagator == "Analytical secular J2"
        @test result_j2.elements.a[1] == result_j2.elements.a[end]
        @test result_j2.elements.e[1] == result_j2.elements.e[end]
        @test result_j2.elements.i[1] == result_j2.elements.i[end]
    end

    @testset "J2 secular rates" begin
        Ωdot, ωdot, Mdot = compute_j2_secular_rates(
            sat.a, sat.e, sat.i
        )
        @test Ωdot > 0  # retrograde orbit at 98 deg
        @test Mdot > 0

        t = 12345.0
        check = element_rate_check(sat, t)
        @test maximum(abs.(check.residual)) < 1e-12
    end

    @testset "Repeat-orbit design" begin
        i = deg2rad(98.0)
        a = design_repeat_orbit(15, 1, i; e = 0.0)
        @test 6500.0 < a < 8000.0
        @test abs(repeat_orbit_error(a, 0.0, i, 15, 1)) < 1e-9
        @test isapprox(repeat_ratio(a, 0.0, i), 15.0; atol = 1e-8)
    end

    @testset "Closure analysis" begin
        i = deg2rad(98.0)
        a = design_repeat_orbit(15, 1, i; e = 0.0)
        repeat_sat = Satellite(
            "Repeat15x1",
            a, 0.0, i, 0.0, 0.0, 0.0
        )
        T = repeat_period(repeat_sat, 15, 1)
        closure = repeat_cycle_closure(
            repeat_sat,
            T;
            n_cycles = 2,
            n_samples = 25,
        )
        @test size(closure.errors) == (2, 25)
        @test length(closure.rms_error) == 2
        @test length(closure.maximum_error) == 2
        @test all(closure.rms_error .>= 0)
    end
end
