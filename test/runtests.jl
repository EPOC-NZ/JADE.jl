using SDDP, Random, JADE, CRC32c, Test

ENV["JADE_DIR"] = joinpath(dirname(@__DIR__), "test")
modeldir = JADE.@JADE_DIR

include("data_validation.jl")

# function checksum(filename, crc = zero(UInt32), blocksize = 16384)
#     open(filename, "r") do f
#         while !eof(f)
#             crc = CRC32c.crc32c(read(f, blocksize), crc)
#         end
#     end
#     return crc
# end

#--------------------------------------------------------------------------
# Deterministic problem for 20 weeks starting from week 45 in 2008,
# scenarios from 2008, 30 iterations
#--------------------------------------------------------------------------
test1 = joinpath(modeldir, "Input", "test1")

# # hash the test1 directory to ensure that the files are the same.
# crc = zero(UInt32)
# for file in sort(readdir(test1))
#     global crc
#     crc = checksum(joinpath(test1, file), crc)
# end
# @test crc == 0xa4e90411

using GLPK
optimizer = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0)

#--------------------------------------------------------------------------
@testset "Building model from input files." begin
    rundata = define_JADE_model("test1")
    rundata.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    JADEmodel = create_JADE_model(rundata, optimizer)
    simulation = define_JADE_simulation("test1")
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb = JADEmodel.sddpm.most_recent_training_results.log[end].bound
    println("JADE LB  = ", lb)
    println("DOASA LB = ", 128312101.499773)
    @test lb ≈ 128312101.499773 atol = 1e4
    #--------------------------------------------------------------------------
    # Simulate from memory
    @testset "Simulating a policy from memory." begin
        results = JADE.simulate(JADEmodel, simulation)
        println(
            "Objective from simulation 1 =  ",
            results[1][1][:stage_objective] + results[1][1][:bellman_term],
        )
        println(
            "Objective from simulation 2 =  ",
            results[2][1][:stage_objective] + results[2][1][:bellman_term],
        )
        @test lb ≈ results[1][1][:stage_objective] + results[1][1][:bellman_term] atol = 1e5
    end

    @testset "Rebuilding model from cuts then simulating." begin
        println("Simulating from cut files.")
        JADEmodel = create_JADE_model(rundata, optimizer)
        results = JADE.simulate(JADEmodel, simulation)
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
    end

    @testset "Doing extra iterations with initial cuts." begin
        optimize_policy!(JADEmodel, solve_options, print_level = 0)
        @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
        println("JADE LB  = ", JADEmodel.sddpm.most_recent_training_results.log[end].bound)
        println("DOASA LB = ", 128312101.499773)
        @test JADEmodel.sddpm.most_recent_training_results.log[end].bound ≈ 128312101.499773 atol =
            1e5
    end
end

#--------------------------------------------------------------------------
# Force load shedding with fake data
@testset "Checking load-shedding equations with fake data." begin
    # Set capacities of DC lines going to/from Haywoods to be 0.
    # The arcs chosen from data are (Hay, SI) and (NI, Hay) so refer to nodes in that order
    rundata = define_JADE_model("test1")
    rundata.use_terminal_mwvs = true
    solve_options = define_JADE_solve_options("test1")
    d = JADEdata(rundata)
    fake_transmission1 = JADE.TransArc(
        0.0,
        0.0,
        Tuple{Float64,Float64}[],#zeros(6, 2),
        Tuple{Float64,Float64}[],#zeros(6, 2),
        0.0,
        0.0,
        d.transmission[(:NI, :HAY)].posoutage,
        d.transmission[(:NI, :HAY)].negoutage,
    )
    fake_transmission2 = JADE.TransArc(
        0.0,
        0.0,
        Tuple{Float64,Float64}[],#zeros(6, 2),
        Tuple{Float64,Float64}[],#zeros(6, 2),
        0.0,
        0.0,
        d.transmission[(:HAY, :SI)].posoutage,
        d.transmission[(:HAY, :SI)].negoutage,
    )
    d.transmission[(:NI, :HAY)] = fake_transmission1
    d.transmission[(:HAY, :SI)] = fake_transmission2
    sddpm = JADE.JADEsddp(d, optimizer)
    JADEmodel = JADE.JADEModel(sddpm, d)
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb = JADEmodel.sddpm.most_recent_training_results.log[end].bound
    println("JADE LB  = ", lb)
    println("DOASA LB = ", 9803982135.425320)
    @test lb ≈ 9803982135.425 atol = 1e5
end

@testset "Testing losses / transmission outages" begin
    rundata = define_JADE_model("test7")
    rundata.use_terminal_mwvs = true
    rundata.scale_objective = 1e6
    rundata.first_week_known = true
    @info("Model defined")
    JADEmodel = create_JADE_model(rundata, optimizer)
    solve_options = define_JADE_solve_options("test7")
    @info("Ready to optimize")
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB = ", lb)
    println("DOASA LB = ", 138197246.8597)
    @test lb ≈ 138197246.8597 atol = 0.001

    rundata = define_JADE_model("test2")
    rundata.use_terminal_mwvs = true
    rundata.scale_objective = 1e6
    rundata.first_week_known = true
    @info("Model defined")
    JADEmodel = create_JADE_model(rundata, optimizer)
    solve_options = define_JADE_solve_options("test2")
    @info("Ready to optimize")
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB = ", lb)
    println("DOASA LB = ", 254678146.80814)
    @test lb ≈ 254678146.80814 atol = 0.001
end

#--------------------------------------------------------------------------
# Scenarios 1997 - 2007, 8 weeks, 400 iterations
# Start from 2008, week 45 so most data can be reused.
@testset "Testing a small stochastic problem." begin
    # Update run data
    rundata = define_JADE_model("test1", run_file = "run2")
    rundata.use_terminal_mwvs = true
    rundata.scale_objective = 1e6
    rundata.first_week_known = true
    solve_options = define_JADE_solve_options("test1", run_file = "run2")
    JADEmodel = create_JADE_model(rundata, optimizer)
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB  = ", lb)
    println("DOASA LB = ", 72_224_064.39)
    @test lb ≈ 72_224_064.39 atol = 500_000
end
#---------------------------------------------------------------------------
# Test DIA implementation
@testset "Testing problem with DIA length 3." begin
    # Update run data
    rundata = define_JADE_model("test1", run_file = "run3")
    rundata.use_terminal_mwvs = true
    rundata.scale_objective = 1e6
    rundata.first_week_known = true
    solve_options = define_JADE_solve_options("test1", run_file = "run3")
    JADEmodel = create_JADE_model(rundata, optimizer)
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB =     ", lb)
    println("DOASA LB was: ", 47_193_083.00)
    @test lb ≈ 47_193_083.00 atol = 500_000
end

@testset "Testing infinite horizon" begin
    # Update run data
    rundata = define_JADE_model("test1", run_file = "run4")
    rundata.scale_objective = 1e6
    rundata.weekly_discounting = false
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    JADEmodel = create_JADE_model(rundata, optimizer)
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 1.26e9 atol = 10_000_000

    JADEmodel = create_JADE_model(rundata, optimizer)
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    solve_options.warmstart_cuts = true
    solve_options.iterations = 1
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 1.26e9 atol = 10_000_000

    # rundata = define_JADE_model("test1", run_file = "run4")
    # rundata.scale_objective = 1e6
    # rundata.discount = 0.0
    # rundata.policy_dir = "finalcuts"
    # rundata.use_terminal_mwvs = false
    # solve_options = define_JADE_solve_options("test1", run_file = "run4")
    # solve_options.warmstart_cuts = true
    # solve_options.finalcutsfromweek = 52
    # solve_options.savedcuts = joinpath("test1", "2008_steady_state")
    # JADEmodel = create_JADE_model(rundata, optimizer)
    # optimize_policy!(JADEmodel, solve_options, print_level = 0)
    # @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    # lb =
    #     JADEmodel.sddpm.most_recent_training_results.log[end].bound *
    #     rundata.scale_objective
    # println("JADE LB =     ", lb)
    # @test lb ≈ 1.26e9 atol = 10_000_000

    rundata = define_JADE_model("test1", run_file = "run4")
    rundata.scale_objective = 1e6
    rundata.weekly_discounting = true
    solve_options = define_JADE_solve_options("test1", run_file = "run4")
    solve_options.iterations = 200
    JADEmodel = create_JADE_model(rundata, optimizer)
    optimize_policy!(JADEmodel, solve_options, print_level = 0)
    @test JADEmodel.sddpm.most_recent_training_results.status == :iteration_limit
    lb =
        JADEmodel.sddpm.most_recent_training_results.log[end].bound *
        rundata.scale_objective
    println("JADE LB =     ", lb)
    @test lb ≈ 9.79e8 atol = 10_000_000
end
