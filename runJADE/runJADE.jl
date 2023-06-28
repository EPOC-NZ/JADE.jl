## This file should be in a directory containing an Input directory, which has a directory
## called <data_dir> containing the JADE input files.

using JADE, JuMP

## Choose your solver
using Gurobi
env = Gurobi.Env()
optimizer = optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)

#using CPLEX
#optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)

#using GLPK
#optimizer = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0)

## Set directory containing the Input / Output subdirectories
ENV["JADE_DIR"] = @__DIR__

## Modify these settings
data_dir = "test1"
run_file = "run"
training = true
simulation = true

## Set up inputs for JADE models
rundata = define_JADE_model(data_dir, run_file = run_file)

## Create JADE model from the runfile data
model = create_JADE_model(rundata, optimizer)

if training
    ## Solve options
    solve_options = define_JADE_solve_options(data_dir, run_file = run_file)

    ## Optimize and store policies
    optimize_policy!(model, solve_options)
end

if simulation
    ## Perform simulation
    sim_settings = define_JADE_simulation(data_dir, run_file = run_file)
    results = simulate(model, sim_settings)
end

## Visualise simulation
JADE.plot_storage(results, data_dir * "-" * runfile * "-" * sim_settings.sim_dir)
JADE.plot_prices(results, data_dir * "-" * runfile * "-" * sim_settings.sim_dir)
