#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

module TestCases

using Test
using JADE

import HiGHS
import JuMP
import SDDP

const MOI = JuMP.MOI

function runtests()
    jade_dir = get(ENV, "JADE_DIR", nothing)
    ENV["JADE_DIR"] = joinpath(dirname(@__DIR__), "test")
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    if jade_dir === nothing
        delete!(ENV, "JADE_DIR")
    else
        ENV["JADE_DIR"] = jade_dir
    end
    return
end

"""
    test_case_1()

Deterministic problem for 20 weeks starting from week 45 in 2008,
scenarios from 2008, 30 iterations
"""
function test_case_1()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1")
    data.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    model = create_JADE_model(data, optimizer)
    simulation = define_JADE_simulation("test1")
    optimize_policy!(model, solve_options; print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound
    DOASA_LB = 128_312_101.499_773
    println("JADE LB  = ", lb)
    println("DOASA LB = ", DOASA_LB)
    @test ≈(lb, DOASA_LB; atol = 1e4)
    # Simulating a policy from memory
    results = JADE.simulate(model, simulation)
    println(
        "Objective from simulation 1 =  ",
        results[1][1][:stage_objective] + results[1][1][:bellman_term],
    )
    println(
        "Objective from simulation 2 =  ",
        results[2][1][:stage_objective] + results[2][1][:bellman_term],
    )
    @test lb ≈ results[1][1][:stage_objective] + results[1][1][:bellman_term] atol = 1e5
    # Rebuilding model from cuts then simulating
    println("Simulating from cut files.")
    model = create_JADE_model(data, optimizer)
    results = JADE.simulate(model, simulation)
    println(
        "Objective from simulation 1 =  ",
        results[1][1][:stage_objective] + results[1][1][:bellman_term],
    )
    println(
        "Objective from simulation 2 =  ",
        results[2][1][:stage_objective] + results[2][1][:bellman_term],
    )
    # 1.284976116274637e8 from new code
    # Something within $2M of LB acceptable
    @test lb ≈ results[1][1][:stage_objective] + results[1][1][:bellman_term] atol = 2e6
    # Doing extra iterations with initial cuts.
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    println("JADE LB  = ", model.sddpm.most_recent_training_results.log[end].bound)
    println("DOASA LB = ", DOASA_LB)
    @test ≈(model.sddpm.most_recent_training_results.log[end].bound, DOASA_LB; atol = 1e5)
    return
end

function test_case_1_load_shedding()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    # Set capacities of DC lines going to/from Haywoods to be 0.
    # The arcs chosen from data are (Hay, SI) and (NI, Hay) so refer to nodes in that order
    data = define_JADE_model("test1")
    data.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    d = JADEdata(data)
    fake_transmission1 = JADE.TransArc(
        0.0,
        0.0,
        Tuple{Float64,Float64}[],
        Tuple{Float64,Float64}[],
        0.0,
        0.0,
        d.transmission[(:NI, :HAY)].posoutage,
        d.transmission[(:NI, :HAY)].negoutage,
    )
    fake_transmission2 = JADE.TransArc(
        0.0,
        0.0,
        Tuple{Float64,Float64}[],
        Tuple{Float64,Float64}[],
        0.0,
        0.0,
        d.transmission[(:HAY, :SI)].posoutage,
        d.transmission[(:HAY, :SI)].negoutage,
    )
    d.transmission[(:NI, :HAY)] = fake_transmission1
    d.transmission[(:HAY, :SI)] = fake_transmission2
    sddpm = JADE.JADEsddp(d, optimizer)
    model = JADE.JADEModel(sddpm, d)
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound
    println("JADE LB  = ", lb)
    DOASA_LB = 9_803_982_135.425_320
    println("DOASA LB = ", DOASA_LB)
    @test ≈(lb, DOASA_LB; atol = 1e5)
    return
end

function test_case_1_stochastic()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1", run_file = "run2")
    data.use_terminal_mwvs = true
    data.scale_objective = 1e6
    data.first_week_known = true
    solve_options = define_JADE_solve_options("test1", run_file = "run2")
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB  = ", lb)
    DOASA_LB = 72_224_064.39
    println("DOASA LB = ", DOASA_LB)
    @test ≈(lb, DOASA_LB; atol = 5e5)
    return
end

function test_case_1_dynamic_inflow_adjustment()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1", run_file = "run3")
    data.use_terminal_mwvs = true
    data.scale_objective = 1e6
    data.first_week_known = true
    solve_options = define_JADE_solve_options("test1", run_file = "run3")
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB =     ", lb)
    DOASA_LB = 47_193_083.00
    println("DOASA LB was: ", DOASA_LB)
    @test ≈(lb, DOASA_LB; atol = 5e5)
    return
end

function test_case_1_infinite_horizon()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    # Update run data
    data = define_JADE_model("test1", run_file = "run4")
    data.scale_objective = 1e6
    data.weekly_discounting = false
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 1.26e9 atol = 10_000_000

    model = create_JADE_model(data, optimizer)
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    solve_options.warmstart_cuts = true
    solve_options.iterations = 1
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 1.26e9 atol = 10_000_000
    data = define_JADE_model("test1", run_file = "run4")
    data.scale_objective = 1e6
    data.weekly_discounting = true
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    solve_options.iterations = 200
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options, print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 9.79e8 atol = 10_000_000
    return
end

function test_case_7()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test7")
    data.use_terminal_mwvs = true
    data.scale_objective = 1e6
    data.first_week_known = true
    @info("Model defined")
    model = create_JADE_model(data, optimizer)
    solve_options = define_JADE_solve_options("test7")
    @info("Ready to optimize")
    optimize_policy!(model, solve_options, print_level = 0)
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB = ", lb)
    println("DOASA LB = ", 138197246.8597)
    @test lb ≈ 138197246.8597 atol = 0.001

    data = define_JADE_model("test2")
    data.use_terminal_mwvs = true
    data.scale_objective = 1e6
    data.first_week_known = true
    @info("Model defined")
    model = create_JADE_model(data, optimizer)
    solve_options = define_JADE_solve_options("test2")
    @info("Ready to optimize")
    optimize_policy!(model, solve_options, print_level = 0)
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB = ", lb)
    DOASA_LB = 254_678_146.808_14
    println("DOASA LB = ", DOASA_LB)
    @test ≈(lb, DOASA_LB; atol = 0.001)
    return
end

function test_case_1_infinite_horizon_then_finite_horizon()
    data = define_JADE_model("test1"; run_file = "run4")
    data.scale_objective = 1e6
    data.weekly_discounting = false
    options = define_JADE_solve_options("test1"; run_file = "run4")
    options.write_eohcuts = true
    model = create_JADE_model(data, HiGHS.Optimizer)
    optimize_policy!(model, options; print_level = 0)
    @test model.sddpm.most_recent_training_results.status == :iteration_limit
    lb = model.sddpm.most_recent_training_results.log[end].bound * data.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 1.26e9 atol = 10_000_000
    # Now build a finite horizon policy using the EOH cuts
    data.steady_state = false
    options.write_eohcuts = false
    options.reset_starting_levels = true
    data.number_of_wks = 10
    options.eoh_cutfile = "2008_steady_state"
    model_finite = create_JADE_model(data, HiGHS.Optimizer)
    num_cons_pre = JuMP.num_constraints(
        model_finite.sddpm[10].subproblem;
        count_variable_in_set_constraints = false,
    )
    optimize_policy!(model_finite, options; print_level = 0)
    @test length(model_finite.sddpm.nodes) == 10
    num_cons_post = JuMP.num_constraints(
        model_finite.sddpm[10].subproblem;
        count_variable_in_set_constraints = false,
    )
    @test num_cons_post > num_cons_pre
    @test model_finite.sddpm.most_recent_training_results.status == :iteration_limit
    @test ≈(
        SDDP.calculate_bound(model_finite.sddpm) * data.scale_objective,
        1.2575e9;
        atol = 1e6,
    )
    # Now build a finite horizon policy using the terminal water value
    options.eoh_cutfile = ""
    data.use_terminal_mwvs = true
    model_finite_no_eoh = create_JADE_model(data, HiGHS.Optimizer)
    num_cons_pre_no_eoh = JuMP.num_constraints(
        model_finite_no_eoh.sddpm[10].subproblem;
        count_variable_in_set_constraints = false,
    )
    optimize_policy!(model_finite_no_eoh, options; print_level = 0)
    @test length(model_finite_no_eoh.sddpm.nodes) == 10
    num_cons_post_no_eoh = JuMP.num_constraints(
        model_finite_no_eoh.sddpm[10].subproblem;
        count_variable_in_set_constraints = false,
    )
    @test num_cons_pre < num_cons_pre_no_eoh  # no_eoh larger because of terminal
    @test num_cons_post > num_cons_post_no_eoh  # eoh larger because of cuts
    @test model_finite_no_eoh.sddpm.most_recent_training_results.status == :iteration_limit
    @test ≈(
        SDDP.calculate_bound(model_finite_no_eoh.sddpm) * data.scale_objective,
        0.3256e9;
        atol = 1e5,
    )
    return
end

function test_case_1_simulate_historical()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1")
    data.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options; print_level = 0)
    simulation = define_JADE_simulation("test1")
    simulation.sim_type = :historical
    @test_throws(
        ErrorException("Invalid settings found. See REPL for details."),
        JADE.simulate(model, simulation),
    )
    simulation.sim_years = [2008]
    @test_throws(
        ErrorException("Invalid settings found. See REPL for details."),
        JADE.simulate(model, simulation),
    )
    simulation.replications = 1
    results = JADE.simulate(model, simulation)
    @test length(results) == 1
    @test length(results[1]) == 20
    @test sum(results[1][1][:lostload]) ≈ 0.0
    @test results[1][1][:inflow][:LAKE_AVIEMORE] == 16.2042285714286
    @test results[1][20][:inflow][:LAKE_AVIEMORE] == 8.64508571428571
    return
end

function test_case_1_simulate_historical_cyclic()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1"; run_file = "run4")
    data.scale_objective = 1e6
    data.weekly_discounting = false
    options = define_JADE_solve_options("test1"; run_file = "run4")
    model = create_JADE_model(data, HiGHS.Optimizer)
    optimize_policy!(model, options; print_level = 0)
    simulation = define_JADE_simulation("test1")
    simulation.sim_type = :historical
    simulation.replications = 2
    simulation.sim_years = [2008, 2009]
    results = JADE.simulate(model, simulation)
    @test length(results) == 2
    @test length(results[1]) == 52
    @test sum(results[1][1][:lostload]) ≈ 0.0
    @test results[1][1][:inflow][:LAKE_AVIEMORE] ≈ 11.9486285714286
    @test results[1][20][:inflow][:LAKE_AVIEMORE] ≈ 8.0728
    @test results[2][1][:inflow][:LAKE_AVIEMORE] ≈ 46.7389142857143
    @test results[2][20][:inflow][:LAKE_AVIEMORE] ≈ 86.2477714285714
    return
end

function test_case_1_simulate_monte_carlo_cyclic()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1"; run_file = "run4")
    data.scale_objective = 1e6
    data.weekly_discounting = false
    options = define_JADE_solve_options("test1"; run_file = "run4")
    model = create_JADE_model(data, HiGHS.Optimizer)
    optimize_policy!(model, options; print_level = 0)
    simulation = define_JADE_simulation("test1")
    simulation.sim_type = :monte_carlo
    simulation.replications = 5
    results = JADE.simulate(model, simulation)
    @test length(results) == 5
    @test length(results[1]) == 52
    @test sum(results[1][1][:lostload]) ≈ 0.0
    # Test that the simulation was cyclic
    res_level = results[1][52][:reslevel][:LAKE_PUKAKI].out
    @test res_level == results[2][1][:reslevel][:LAKE_PUKAKI].in
    return
end

function test_case_1_simulate_historical_initial_state()
    optimizer = MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true)
    data = define_JADE_model("test1")
    data.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    model = create_JADE_model(data, optimizer)
    optimize_policy!(model, solve_options; print_level = 0)
    simulation = define_JADE_simulation("test1")
    simulation.sim_type = :historical
    simulation.replications = 1
    simulation.sim_years = [2008]
    results = JADE.simulate(model, simulation)
    # Test that setting the initial state as the default gives the same answer
    simulation.initial_state =
        Dict(String(k) => v for (k, v) in model.sddpm.initial_root_state)
    results2 = JADE.simulate(model, simulation)
    @test results[1][1][:total_storage] == results2[1][1][:total_storage]
    delete!(simulation.initial_state, "reslevel[LAKE_OHAU]")
    @test_throws(
        ErrorException("initial_state dictionary has incorrect number of states"),
        JADE.simulate(model, simulation),
    )
    simulation.initial_state["reslevel[LAKE_FAKE]"] = 0.0
    @test_throws(
        ErrorException("Invalid state: reslevel[LAKE_FAKE]"),
        JADE.simulate(model, simulation),
    )
    return
end

end  # module

TestCases.runtests()
