#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

struct JADEModel
    sddpm::SDDP.PolicyGraph
    d::JADEData
end

function read_string(str::String)
    if length(str) <= 1
        return str
    elseif string(str[1]) == "\"" && string(str[end]) == "\""
        return str[2:end-1]
    else
        return str
    end
end

"""
    define_JADE_model(inputdir::String; run_file::String = "run")

This function reads a JADE run file in the `Input/<inputdir>` folder. In particular,
this function reads the parameters that are needed in order to define the SDDP model.

### Required Arguments
`inputdir` is the name of the subdirectory (within `Input`) that contains the JADE
data files.

### Keyword Arguments
`run_file` is the name of the csv file that contains the parameters we wish to load.
"""
function define_JADE_model(inputdir::String; run_file::String = "run")
    file = joinpath(@__JADE_DIR__, "Input", inputdir, run_file * ".csv")
    run_options = parse_run_options(file)

    # Scaling factor for reservoir volume
    scale_res = 1.0
    if haskey(run_options, "Scale reservoirs")
        scale_res = parse(Float64, run_options["Scale reservoirs"])
    end

    scale_obj = 1.0
    if haskey(run_options, "Scale objective")
        scale_obj = parse(Float64, run_options["Scale objective"])
    end

    ss = 0.0
    if haskey(run_options, "Steady state")
        ss = parse(Float64, run_options["Steady state"])
    end

    if haskey(run_options, "System")
        @warn("'System' keyword is used by DOASA, not JADE. It has been ignored.")
    end

    losses = :default
    if haskey(run_options, "Loss model")
        if replace(run_options["Loss model"], "\"" => "") ∈
           ["piecewise linear", "piecewise"]
            losses = :piecewise
        elseif replace(run_options["Loss model"], "\"" => "") == "none"
            losses = :none
        elseif replace(run_options["Loss model"], "\"" => "") == "resistive"
            losses = :resistive
        else
            error("Invalid loss model: " * run_options["Loss model"])
        end
    end

    if haskey(run_options, "Policy name") && run_options["Policy name"] != "\"\""
        policy_dir = replace(run_options["Policy name"], "\"" => "")
    else
        policy_dir = inputdir
    end

    sample_start = parse(Int, run_options["Sample start year"])
    sample_end = parse(Int, run_options["Sample end year"])

    scenario_dir =
        haskey(run_options, "Scenario") ? read_string(run_options["Scenario"]) : ""

    return define_JADE_model(
        data_dir = inputdir,
        policy_dir = policy_dir,
        start_year = parse(Int, run_options["Problem start year"]),
        start_week = parse(Int, run_options["Problem start week"]),
        total_weeks = parse(Int, run_options["Number of weeks"]),
        hydro_sample_range = collect(sample_start:sample_end),
        discount = ss,
        losses = losses,
        DIA = parse(Int, run_options["Inflow correlation length"]),
        flow_penalties = (
            parse(Float64, run_options["LB flow penalty"]),
            parse(Float64, run_options["UB flow penalty"]),
        ),
        scale_reservoirs = scale_res,
        scale_objective = scale_obj,
        scenario_dir = scenario_dir,
    )
end

"""
    define_JADE_model(;
        data_dir::String = "",
        policy_dir::String = "",
        start_year::Int = 0,
        start_week::Int = 1,
        total_weeks::Int = 52,
        hydro_sample_range::Union{Tuple{Int,Int},Int,UnitRange{Int},Vector} = 0,
        discount::Real = 0.0,
        losses::Symbol = :default,
        DIA::Int = 1,
        flow_penalties::Union{Tuple{Real,Real},NamedTuple{(:under, :over),<:Tuple{Real,Real}}} = (
            500.0,
            50.0,
        ),
        weekly_discounting::Bool = true,
        scale_reservoirs::Real = 1.0,
        scale_objective::Real = 1.0,
        use_terminal_mwvs::Bool = false,
        first_week_known::Bool = false,
        steady_state::Union{Nothing,Bool} = nothing,
        decision_rules::Vector{DecisionRule} = DecisionRule[],
        scenario_dir::String = "",
    )

This function sets the parameters that are needed in order to define the SDDP model.

### Required Arguments
`data_dir` is the name of the subdirectory (within `Input`) that contains the JADE
data files.

`policy_dir` is the name of the subdirectory (within `Output/data_dir`) that the cuts
and simulation output will be stored in.

`start_year` is the year that the policy is created for.

`hydro_sample_range` is the set of years that define the hydro inflow distribution.

### Keyword Arguments
`start_week` is the first week of the `start_year`.

`total_weeks` is the total number of weeks that the policy is created for.

`discount` is the discount factor used for the steady-state model, must be less than 1.
Setting this to 0 will set JADE to a finite horizon model. The terminal water values will
only be used if `discount` is 0.

`losses` sets the loss model; this can be `:none`, `:piecewise` or `:resistive`.

`DIA` is the dynamic inflow adjustment parameter; 0 or 1 means no adjustment, higher values
simulate an increased temporal correlation by increasing the variance.

`flow_penalties` is a Tuple of the cost of violating the lower bound, upper bound flow constraints,
respectively.

`weekly_discounting` is a `Bool` which sets the discounting to occur in each stage rather than at the end of the horizon.

`scale_reservoirs` this may be increased to improve numerical stability (~100).

`scale_objective` scales down the objective function, this also scales down cuts (~1000).

`use_terminal_mwvs` if set to `true` this will use the MWVs stored in `terminal_water_value.csv`, if set to `false`
these won't be used.

`first_week_known` if set to `true` the training will use deterministic inflows for the first week of the model. These will be
sourced from the inflows recorded for that week and year.

`steady_state` if set to `true` the training will operate in steady-state mode, if `false` the model will use a finite-horizon.

`decision_rules` if any reservoir has restrictions on output, as a function of the storage, a `DecisionRule` can be applied. A vector
of such decisions can be supplied via this argument.
"""
function define_JADE_model(;
    data_dir::String = "",
    policy_dir::String = "",
    start_year::Int = 0,
    start_week::Int = 1,
    total_weeks::Int = 52,
    hydro_sample_range::Union{Tuple{Int,Int},Int,UnitRange{Int},Vector} = 0,
    discount::Real = 0.0,
    losses::Symbol = :default,
    DIA::Int = 1,
    flow_penalties::Union{Tuple{Real,Real},NamedTuple{(:under, :over),<:Tuple{Real,Real}}} = (
        500.0,
        50.0,
    ),
    weekly_discounting::Bool = true,
    scale_reservoirs::Real = 1.0,
    scale_objective::Real = 1.0,
    use_terminal_mwvs::Bool = false,
    first_week_known::Bool = false,
    steady_state::Union{Nothing,Bool} = nothing,
    decision_rules::Vector{DecisionRule} = DecisionRule[],
    scenario_dir::String = "",
)
    errors = ""

    if typeof(flow_penalties) <: Union{NamedTuple,Tuple}
        if typeof(flow_penalties) <: NamedTuple
            flow_penalties = (flow_penalties[:under], flow_penalties[:over])
        end
    else
        errors *= "\n'flow_penalties' should be a NamedTuple: (under = x, over = y)."
    end

    if errors != ""
        @error("Error defining JADE model:" * errors)
        error("Error defining JADE model. See REPL for details.")
    end

    if steady_state == nothing
        if discount == 0.0
            steady_state = false
        else
            steady_state = true
        end
    end

    hydro_sample_range = convert_to_int_array(hydro_sample_range)

    rundata = RunData(
        data_dir,
        start_year,
        start_week,
        total_weeks,
        hydro_sample_range,
        DIA,
        flow_penalties[1],
        flow_penalties[2],
        policy_dir,
        discount,
        losses,
        length(hydro_sample_range),
        weekly_discounting,
        scale_reservoirs,
        scale_objective,
        use_terminal_mwvs,
        first_week_known,
        steady_state,
        decision_rules,
        scenario_dir,
    )

    check_settings_compatibility(rundata = rundata)

    return rundata
end

"""
    define_JADE_simulation(inputdir::String; run_file = "run")

This function reads a csv file containing the settings for a JADE simulation.

### Required Arguments
`inputdir` is the name of the subdirectory (within `Input`) where the `run_file`
is located.

### Keyword Arguments
`run_file` is the name of the csv file containing the parameters for the simulation.
"""
function define_JADE_simulation(inputdir::String; run_file = "run")
    file = joinpath(@__JADE_DIR__, "Input", inputdir, run_file * ".csv")
    run_options = parse_run_options(file)

    reset_starting_levels = :default

    if haskey(run_options, "Reset initial simulation state")
        reset_starting_levels =
            parse(Bool, lowercase(run_options["Reset initial simulation state"]))
    end

    JADEsim = nothing

    if haskey(run_options, "Simulation name")
        simulation_name = read_string(run_options["Simulation name"])
    else
        @warn("No 'Simulation name' specified, defaulting to $(inputdir)")
        simulation_name = inputdir
    end

    if lowercase(read_string(run_options["Simulation type"])) == "monte carlo"
        JADEsim = define_JADE_simulation(
            sim_dir = simulation_name,
            sim_type = :monte_carlo,
            replications = parse(Int, run_options["Simulation sample size"]),
            reset_starting_levels = reset_starting_levels,
            random_seed = parse(Int, run_options["Random seed"]),
        )
    elseif lowercase(read_string(run_options["Simulation type"])) == "historical"
        sample_start =
            haskey(run_options, "Simulation start year") ?
            parse(Int, run_options["Simulation start year"]) :
            parse(Int, run_options["Sample start year"])

        sample_end =
            haskey(run_options, "Simulation end year") ?
            parse(Int, run_options["Simulation end year"]) :
            parse(Int, run_options["Sample end year"])

        sample_size = 0
        randomize_years = false
        if haskey(run_options, "Simulation sample size")
            randomize_years = true
            sample_size = parse(Int, run_options["Simulation sample size"])
        end

        random_seed =
            haskey(run_options, "Random seed") ? parse(Int, run_options["Random seed"]) :
            Int(floor((time() - floor(time())) * 10000))

        JADEsim = define_JADE_simulation(
            sim_dir = simulation_name,
            sim_type = :historical,
            replications = sample_size,
            sim_years = sample_start:sample_end,
            randomize_years = randomize_years,
            reset_starting_levels = reset_starting_levels,
            random_seed = random_seed,
        )
    else
        @error("Undefined simulation type")
    end

    return JADEsim
end

function define_JADE_simulation(
    sim_dir::String,
    sim_type::Symbol,
    replications::Int = 0;
    sim_years = nothing,
    randomize_years::Union{Symbol,Bool} = :default,
    reset_starting_levels::Union{Symbol,Bool} = :default,
    number_of_cycles::Int = 1,
    initial_stage::Int = 1,
    initial_state::Union{Nothing,Dict{String,<:Real}} = nothing,
    random_seed::Int = Int(floor((time() - floor(time())) * 10000)),
)
    @warn(
        "This version of define_JADE_simulation function is only provided for backward compatibility."
    )

    return define_JADE_simulation(;
        sim_dir = sim_dir,
        sim_type = sim_type,
        replications = replications,
        sim_years = sim_years,
        randomize_years = randomize_years,
        reset_starting_levels = reset_starting_levels,
        number_of_cycles = number_of_cycles,
        initial_stage = initial_stage,
        initial_state = initial_state,
        random_seed = random_seed,
    )
end

"""
    define_JADE_simulation(;
        sim_dir::String = "",
        sim_type::Symbol = :not_set,
        replications::Int = 0,
        sim_years = nothing,
        randomize_years::Bool = false,
        reset_starting_levels::Bool = false,
        number_of_cycles::Int = 1,
        initial_stage::Int = 1,
        initial_state::Union{Nothing,Dict{String,<:Real}} = nothing,
        random_seed::Int = Int(floor((time() - floor(time())) * 10000)),
    )

This function defines the settings for a JADE simulation.

### Keyword Arguments
`sim_dir` is the name of the subdirectory (within `Output/<data_dir>/<policy_dir>`)
in which the simulation output will be stored.

`sim_type` is the type of simulation; either `:monte_carlo` or `:historical`.

`replications` is the number of sequences of inflows that are simulated. If the
`sim_type` is `:monte_carlo` or `randomize_years` set to `true` then this
argument must be set, otherwise it must not be set.

`sim_years` specifies the years that will be simulated.

`randomize_years` is used with historical simulations; if `true` this causes the
simulation to randomly sample from `sim_years`, whereas if `false` the simulation
will sample the years in sequence.

`reset_starting_levels` is used to reset the reservoirs in a steady-state model.
This will enable a simulation to sample a number of inflows sequences, with the
reservoirs starting at the from the same points for each sequence.

`number_of_cycles` is used with `reset_starting_levels` to show a simulation of multiple
years, from a single single starting point. Must be used with a steady-state
policy.

`initial_state` is a vector of initial storage levels for the simulation, indexed by the reservoirs.

`initial_stage` is the stage of the trained model that the simulation will start from.

`random_seed` sets the random seed prior to generating the random inflow sequences.
"""
function define_JADE_simulation(;
    sim_dir::String = "",
    sim_type::Symbol = :not_set,
    replications::Int = 0,
    sim_years = nothing,
    randomize_years::Union{Symbol,Bool} = :default,
    reset_starting_levels::Union{Symbol,Bool} = :default,
    number_of_cycles::Int = 1,
    initial_stage::Int = 1,
    initial_state::Union{Nothing,Dict{String,<:Real}} = nothing,
    random_seed::Int = Int(floor((time() - floor(time())) * 10000)),
)
    if reset_starting_levels == :default
        if number_of_cycles > 1
            reset_starting_levels = true
        end
    end

    if initial_state != nothing
        temp = Dict{String,Float64}()
        for (k, s) in initial_state
            temp["reslevel["*k*"]"] = s
        end
        initial_state = temp
    elseif initial_stage != 1
        @warn("initial_stage has been set, but initial_state has not")
    end

    if sim_type == :monte_carlo && randomize_years != :default
        @warn(
            "You cannot use 'randomize_years' with :monte_carlo simulations. This will be ignored."
        )
        randomize_years = false
    elseif randomize_years == :default
        randomize_years = false
    end

    if replications != 0 && (sim_type == :historical && !randomize_years)
        @warn(
            "You should not specify the number of 'replications' for sequential historical simulations. The input value has been overwritten."
        )
    end

    if sim_type == :historical
        if sim_years != nothing
            sim_years = convert_to_int_array(sim_years)
            if randomize_years
                Random.seed!(random_seed)

                years = Int[]
                for i in 1:replications*number_of_cycles
                    push!(years, sim_years[Int(floor(rand() * length(sim_years)))+1])
                end

                sim_years = years
            else
                replications = length(sim_years)
                if number_of_cycles > 1
                    years = Int[]
                    for y in sim_years
                        for i in 0:(number_of_cycles-1)
                            push!(years, y + i)
                        end
                    end
                    sim_years = years
                end
            end
        end
    end

    simulation = JADESimulation(
        sim_dir,
        sim_type,
        replications,
        sim_years,
        randomize_years,
        reset_starting_levels,
        number_of_cycles,
        initial_stage,
        initial_state,
        random_seed,
    )

    check_settings_compatibility(simulation = simulation)

    return simulation
end

"""
    create_JADE_model(rundata::RunData, optimizer)

This function creates a `JADEModel` consisting of an SDDP policy graph, and the
associated data.

### Required Arguments
`rundata` is the `RunData` object containing all the SDDP model details.

`optimizer` is the JuMP optimizer you wish to use to train the policy.
"""
function create_JADE_model(rundata::RunData, optimizer = nothing; solver = nothing)
    if solver != nothing
        @warn(
            "The 'solver =' option is provided only for backward compatibility.\n Keyword arguments are no longer needed for the optimizer."
        )
        optimizer = solver
    end

    check_settings_compatibility(rundata = rundata)

    @info("Loading data files...")
    d = JADEdata(rundata)

    @info("Building SDDP model...")
    sddpm = JADEsddp(d, optimizer)

    @info("JADE model created")
    return JADEModel(sddpm, d)
end

"""
    define_JADE_solve_options(inputdir::String; run_file = "run")

This function loads the training options for JADE.

### Required Arguments
`inputdir` is the name of the subdirectory (within `Input`) where the `run_file`
is located.

### Keyword Arguments
`run_file` is the name of the csv file containing the parameters for the simulation.
"""
function define_JADE_solve_options(inputdir::String; run_file = "run")
    file = joinpath(@__JADE_DIR__, "Input", inputdir, run_file * ".csv")
    run_options = parse_run_options(file)

    riskmeasure = (0.0, 1.0)
    if haskey(run_options, "Risk lambda") && haskey(run_options, "Risk beta")
        β = parse(Float64, run_options["Risk beta"])
        λ = parse(Float64, run_options["Risk lambda"])
        if abs(β - 1) > 1e-7 && abs(λ - 1) > 1e-7
            riskmeasure = (λ, β)
        end
    end

    ss = 0.0
    if "Steady state" ∈ keys(run_options)
        ss = parse(Float64, run_options["Steady state"])
    end

    warm = false
    saved_cuts = ""
    if "Use saved cuts from" in keys(run_options)
        saved_cuts = run_options["Use saved cuts from"]
        saved_cuts = replace(saved_cuts, "\"" => "")
        warm = saved_cuts == "" ? false : true
    end

    cutselection = 1
    if "Cut selection" in keys(run_options)
        cutselection = parse(Int, run_options["Cut selection"])
    end

    fractionMC = 1.0
    if "Fraction of forward passes using Monte Carlo inflow sequences" in keys(run_options)
        fractionMC = parse(
            Float64,
            run_options["Fraction of forward passes using Monte Carlo inflow sequences"],
        )
    end

    custom_inflow_file = ""
    if "Custom inflow sequences from" in keys(run_options)
        custom_inflow_file = run_options["Custom inflow sequences from"]
        custom_inflow_file = replace(custom_inflow_file, "\"" => "")
    end

    eoh_cutfile =
        haskey(run_options, "Load EOH cuts file") ?
        replace(run_options["Load EOH cuts file"], "\"" => "") : ""
    write_eohcuts =
        haskey(run_options, "Write EOH cuts file") ?
        parse(Bool, lowercase(run_options["Write EOH cuts file"])) : false

    random_seed =
        haskey(run_options, "Random seed") ? parse(Int, run_options["Random seed"]) :
        Int(floor((time() - floor(time())) * 10000))

    return define_JADE_solve_options(
        iterations = parse(Int, run_options["Maximum iterations"]),
        random_seed = random_seed,
        risk = riskmeasure,
        warmstart_cuts = warm,
        saved_cuts = saved_cuts,
        cutselection = cutselection,
        fractionMC = fractionMC,
        custom_inflow_file = custom_inflow_file,
        eoh_cutfile = eoh_cutfile,
        write_eohcuts = write_eohcuts,
    )
end

"""
    define_JADE_solve_options(;
        iterations::Int = 0,
        random_seed::Int = Int(floor((time() - floor(time())) * 10000)),
        risk::Union{Tuple{Real,Real},NamedTuple{(:lambda, :beta),<:Tuple{Real,Real}}} = (
            0.0,
            1.0,
        ),
        warmstart_cuts::Bool = false,
        reset_starting_levels::Union{Symbol,Bool} = :default,
        saved_cuts::String = "",
        cutselection::Int = 1,
        fractionMC::Float64 = 1.0
        custom_inflow_file::String = ""
        loadcuts_if_nonempty::Bool = false,
        eoh_cutfile::String = "",
        write_eohcuts::Bool = false
    )

This function defines the training options for JADE.

### Keyword Arguments

`iterations` is the number of iterations of SDDP to complete when `optimize_policy!`
is called.

`random_seed` sets the random seed prior to sampling for the forward pass.

`risk` is a Tuple `(λ,β)` defining the objective weight `λ` and the tail
probability `β`.

`warmstart_cuts` if set to `true`, this will load stored cuts before training.

`reset_starting_levels` is used to reset the reservoirs in a steady-state model.
This will enable a steady-state model to focus on a single starting set of
reservoir levels (after initially finding an end-of-horizon value function).

`saved_cuts` is the policy subdirectory from which the cuts should be loaded. If
not specified, cuts will be loaded from the policy subdirectory set for the model.

`cutselection` is an SDDP.jl parameter.

`fractionMC` is the proportion of forward passes that are sampled using Monte-Carlo;
the rest are sequentially sampled from the `sequences`.

`custom_inflow_file` is the name of the file in the `Input/data_dir` directory containing
the custom inflow sequences.

`loadcuts_if_nonempty` if set to `true` cuts will be loaded from file even if cuts are
already in the model.

`eoh_cutfile` the name of end-of-horizon cut file for finite-horizon models.

`write_eohcuts` if set to `true` steady-state models will write end-of-horizon cuts files.

"""
function define_JADE_solve_options(;
    iterations::Int = 0,
    random_seed::Int = Int(floor((time() - floor(time())) * 10000)),
    risk::Union{Tuple{Real,Real},NamedTuple{(:lambda, :beta),<:Tuple{Real,Real}}} = (
        0.0,
        1.0,
    ),
    warmstart_cuts::Bool = false,
    reset_starting_levels::Union{Symbol,Bool} = :default,
    saved_cuts::String = "",
    cutselection::Int = 0,
    fractionMC::Real = 1.0,
    custom_inflow_file::String = "",
    loadcuts_if_nonempty::Bool = false,
    eoh_cutfile::String = "",
    write_eohcuts::Bool = false,
)
    errors = ""

    if typeof(risk) <: Union{NamedTuple,Tuple}
        if typeof(risk) <: NamedTuple
            risk = (risk[:lambda], risk[:beta])
        end
    else
        errors *= "'risk' should be a NamedTuple: (lambda = x, beta = y).\n"
    end

    if risk[1] < 0.0 || risk[1] > 1.0
        errors *= "λ (risk[:lambda] / risk[1]) must be ⊂ [0,1]\n"
    elseif risk[1] == 1.0 && risk[2] < 1.0
        @warn(
            "λ (risk[:lambda] / risk[1]) is 1.0; no weight has been placed on expectation."
        )
    end

    if risk[1] > 0.0 && (risk[2] <= 0.0 || risk[2] > 1.0)
        errors *= "β (risk[:beta] / risk[2]) must be ⊂ (0,1]\n"
    end

    if warmstart_cuts && eoh_cutfile != ""
        @warn(
            "The 'eoh_cutfile': $(eoh_cutfile) will not be used because cuts are being warm_started."
        )
    end

    if !warmstart_cuts && saved_cuts != ""
        @warn(
            "The 'saved_cuts': $(saved_cuts) will not be used because cuts are not being warm_started."
        )
    end

    if errors != ""
        error(errors)
    end

    if risk[1] > 0.0 && risk[2] < 1.0
        # SDDP.jl λ is backwards and puts weight on expected value
        riskmeasure = SDDP.EAVaR(
            lambda = (1 - convert(Float64, risk[1])),
            beta = convert(Float64, risk[2]),
        )
    else
        riskmeasure = SDDP.Expectation()
    end

    if custom_inflow_file == ""
        if fractionMC < 1.0
            custom_inflow_file = "custom_sequences.csv"
        end
    else
        if fractionMC == 1.0
            @warn(
                "Custom inflow file specified, but fractionMC is 1.0. Custom inflows will be ignored."
            )
        end
    end

    solveoptions = JADESolveOptions(
        iterations,
        riskmeasure,
        warmstart_cuts,
        reset_starting_levels,
        random_seed,
        saved_cuts,
        cutselection,
        fractionMC,
        custom_inflow_file,
        loadcuts_if_nonempty,
        eoh_cutfile,
        write_eohcuts,
    )

    check_settings_compatibility(solveoptions = solveoptions)

    return solveoptions
end

function load_model_parameters(model::String, policy::String)
    return load_model_parameters(joinpath(model, policy))
end

function load_model_parameters(path::String)
    if !isfile(path)
        path = joinpath(@__JADE_DIR__, "Output", path, "rundata.json")
    end

    if !ispath(path)
        error("Parameters for JADE model not found: $path.")
    end
    data = JSON.parsefile(path)

    decisionrules = DecisionRule[]

    for i in 1:length(data["decision_rules"])
        push!(
            decisionrules,
            JADE.DecisionRule(
                data["decision_rules"][i]["slope"],
                data["decision_rules"][i]["intercept"],
                Symbol(data["decision_rules"][i]["boundtype"]),
                Symbol(data["decision_rules"][i]["station"]),
                Symbol(data["decision_rules"][i]["reservoir"]),
                Symbol(data["decision_rules"][i]["flowtype"]),
                data["decision_rules"][i]["weeks"],
            ),
        )
    end

    #return data
    return RunData(
        data["data_dir"],
        data["start_yr"],
        data["start_wk"],
        data["number_of_wks"],
        data["sample_years"],
        data["dialength"],
        data["penalty_lb"],
        data["penalty_ub"],
        data["policy_dir"],
        data["discount"],
        Symbol(data["losses"]),
        data["nscenarios"],
        data["weekly_discounting"],
        data["scale_reservoirs"],
        data["scale_objective"],
        data["use_terminal_mwvs"],
        data["first_week_known"],
        data["steady_state"],
        decisionrules,
        data["scenario_dir"],
    )
end

function load_solve_parameters(model::String, policy::String)
    path = joinpath(@__JADE_DIR__, "Output", model, policy, "solveoptions.json")

    if !ispath(path)
        error("Solve options not found in Output" * joinpath("Output", model, policy))
    end
    data = JSON.parsefile(path)
    risk = data["riskmeasure"]
    if risk[1] > 0.0 && risk[2] < 1.0
        riskmeasure = SDDP.EAVaR(lambda = (1 - risk[1]), beta = risk[2])
    else
        riskmeasure = SDDP.Expectation()
    end

    return JADESolveOptions(
        data["iterations"],
        riskmeasure,
        data["warmstart_cuts"],
        data["reset_starting_levels"],
        data["seed"],
        data["savedcuts"],
        data["cutselection"],
        data["fractionMC"],
        data["custom_inflow_file"],
        data["loadcuts_if_nonempty"],
        data["eoh_cutfile"],
        data["write_eohcuts"],
    )
end

function load_simulation_parameters(model::String, policy::String, simulation::String)
    path =
        joinpath(@__JADE_DIR__, "Output", model, policy, simulation, "sim_parameters.json")
    if !ispath(path)# || !ispath(joinpath(data_dir,rundata.json))
        error("Parameters for simulation not found")
    end
    data = JSON.parsefile(path)
    if haskey(data, "sim_years") &&
       data["sim_years"] != nothing &&
       length(data["sim_years"]) > 0
        data["sim_years"] = Int.(data["sim_years"])
    end
    if typeof(data["reset_starting_levels"]) != Bool
        data["reset_starting_levels"] = Symbol(data["reset_starting_levels"])
    end

    return JADESimulation(
        data["sim_dir"],
        Symbol(data["sim_type"]),
        data["replications"],
        data["sim_years"],
        data["randomize_years"],
        data["reset_starting_levels"],
        data["number_of_cycles"],
        data["initial_stage"],
        data["initial_state"],
        data["random_seed"],
    )
end

"""
    function convert_to_int_array(
        values::T where {T<:Union{Tuple{Int,Int},Int,UnitRange{Int},Vector}},
    )

This function converts various ways of defining a set of integers into a vector of
integers.
"""
function convert_to_int_array(
    values::T where {T<:Union{Tuple{Int,Int},Int,UnitRange{Int},Vector}},
)
    if typeof(values) == Int
        return [values]
    elseif typeof(values) == Tuple{Int,Int}
        if values[2] < values[1]
            error("First value of Tuple cannot be greater than second value")
        end
        return collect(values[1]:values[2])
    elseif typeof(values) == UnitRange{Int}
        temp = collect(values)
        if length(temp) > 0
            return collect(values)
        end
    elseif typeof(values) == Array{Int,1}
        return values
    elseif typeof(values) <: Array
        out = Int[]
        for i in values
            if typeof(i) == Int
                push!(out, i)
            elseif typeof(i) == UnitRange{Int}
                for j in i
                    push!(out, j)
                end
            else
                error("Invalid data type - elements must be Int or UnitRange{Int}")
            end
        end
        return out
    end
end

"""
    check_rundata(rd1::RunData, rd2::RunData, match::Symbol)

This function compares the settings within two different RunData objects. If
`match` is `:full` then JADE will throw an error for any difference, if `:match` is `:EOH` it
will only throw an error for major differences (scaling), otherwise an error will also be throw
if `number_of_wks` is different.
"""
function check_rundata(rd1::RunData, rd2::RunData, match::Symbol)
    if match == :full
        if rd1.start_yr != rd2.start_yr || rd1.start_wk != rd2.start_wk
            error(
                "The starting week / year for the current model do not match the previous model.",
            )
        elseif rd1.sample_years != rd2.sample_years
            error(
                "The training years for the current model do not match the previous model.",
            )
        elseif rd1.dialength != rd2.dialength
            error("The dialength for the current model does not match the previous model.")
        elseif rd1.penalty_lb != rd2.penalty_lb || rd1.penalty_ub != rd2.penalty_ub
            error(
                "The flow penalties for the current model do not match the previous model.",
            )
        elseif rd1.discount != rd2.discount
            error(
                "The discount factor for the current model does not match the previous model.",
            )
        elseif rd1.losses != rd2.losses
            error("The losses for the current model does not match the previous model.")
        elseif rd1.weekly_discounting != rd2.weekly_discounting
            error(
                "The discounting mode for the current model does not match the previous model.",
            )
        elseif rd1.use_terminal_mwvs != rd2.use_terminal_mwvs
            error(
                "The use of terminal MWVs for the current model does not match the previous model.",
            )
        elseif rd1.first_week_known != rd2.first_week_known
            error(
                "The first week's inflows for the current model do not match the previous model.",
            )
        elseif rd1.decision_rules != rd2.decision_rules
            error(
                "The decision rules for the current model do not match the previous model",
            )
        end
    end

    if rd1.scale_reservoirs != rd2.scale_reservoirs ||
       rd1.scale_objective != rd2.scale_objective
        error("The scaling factors for the current model do not match the previous model.")
    elseif match != :eoh && rd1.number_of_wks != rd2.number_of_wks
        error(
            "The number of weeks for the current model does not match the previous model.",
        )
    end
end

function duplicate_data(data::TimeSeries{Dict{Symbol,Float64}}, set)
    temp = Dict{Tuple{Symbol,Symbol},Float64}[]
    for k in keys(data)
        temp2 = Dict{Tuple{Symbol,Symbol},Float64}()
        for s in keys(data[k])
            for b in set
                temp2[(s, b)] = data[k][s]
            end
        end
        push!(temp, temp2)
    end
    return JADE.TimeSeries(data.startpoint, temp)
end

"""
    check_settings_compatibility(;
        rundata::Union{Nothing,RunData} = nothing,
        solveoptions::Union{Nothing,JADESolveOptions} = nothing,
        simulation::Union{Nothing,JADESimulation} = nothing,
    )

This function checks the JADE settings for compatibility / consistency.

### Keyword Arguments

`rundata` contains the JADE model options.
`solveoptions` contains the JADE training options.
`simulation` contains the JADE simulation settings.

"""
function check_settings_compatibility(;
    rundata::Union{Nothing,RunData} = nothing,
    solveoptions::Union{Nothing,JADESolveOptions} = nothing,
    simulation::Union{Nothing,JADESimulation} = nothing,
)
    errors = ""

    if rundata != nothing
        if rundata.data_dir == ""
            errors *= "\nNo 'data_dir' specified. This should be in the Input directory."
        end

        if rundata.policy_dir == ""
            errors *= "\nNo 'policy_dir' specified. This directory will store the cuts for the policy."
        end

        if rundata.start_yr == 0
            errors *= "\nNo 'start_year' specified, this should be the year (in terms of demand / supply) that the policy is built for."
        end

        if rundata.start_wk > 52
            errors *= "\n'start_week' must be less than or equal to 52."
        elseif rundata.start_wk < 1
            errors *= "\n'start_week' must be greater than or equal to 1."
        end

        if rundata.sample_years == [0]
            errors *= "\n'hydro_sample_range' must be specified as an Int, Tuple, Array or a UnitRange."
        end

        if rundata.losses ∉ [:default, :none, :piecewise, :resistive]
            errors *= "\nInvalid loss model; 'loss' must be :none, :piecewise or :resistive."
        end

        if rundata.discount < 0.0 || rundata.discount >= 1.0
            errors *= "\n'discount' must be ∈ [0.0,1.0)."
        end

        if rundata.dialength < 0
            errors *= "\n'DIA' must be >= 0."
        end

        if rundata.discount == 0.0 && rundata.steady_state
            errors *= "\nCannot create a steady-state model with no discounting."
        end

        if rundata.steady_state
            if rundata.number_of_wks != 52
                errors *= "\nSteady-state models can only be run over periods of 52 weeks."
            end
            if rundata.first_week_known
                errors *= "\nThe first week inflows must be stochastic for steady-state models."
            end
        end

        if rundata.use_terminal_mwvs && rundata.steady_state
            errors *= "\nYou cannot use terminal MWVs in steady-state models."
        end
    end

    if solveoptions != nothing
        if solveoptions.iterations == 0
            errors *= "\n'iterations' must be specified."
        elseif solveoptions.iterations < 0
            errors *= "\n'iterations' must be a positive integer."
        end

        if solveoptions.fractionMC < 0.0 || solveoptions.fractionMC > 1.0
            errors *= "\nfractionMC must be between 0 and 1."
        end
    end

    if rundata != nothing && solveoptions != nothing
        if !rundata.steady_state &&
           !rundata.use_terminal_mwvs &&
           solveoptions.eoh_cutfile == "" &&
           solveoptions.savedcuts == ""
            errors *= "\nA finite-horizon model must either have 'use_terminal_mwvs' set to 'true' or have an 'eoh_cutfile' specified."
        end

        if !rundata.steady_state && solveoptions.write_eohcuts
            errors *= "\nA finite-horizon model cannot write an end-of-horizon cuts file."
        end

        if !rundata.steady_state && solveoptions.reset_starting_levels == false
            errors *= "\nA finite-horizon model must reset its starting storage levels each iteration."
        end

        if solveoptions.reset_starting_levels == :default
            if rundata.steady_state
                solveoptions.reset_starting_levels = false
            else
                solveoptions.reset_starting_levels = true
            end
        end

        if rundata.steady_state && solveoptions.eoh_cutfile != ""
            errors *= "\nA steady-state model cannot use an 'eoh_cutfile'."
        end
    end

    if simulation != nothing
        if simulation.sim_dir == ""
            errors *= "\nYou must specify a 'sim_dir', where the simulation output will be stored."
        end

        if simulation.sim_type != :monte_carlo && simulation.sim_type != :historical
            errors *= "\n'sim_type' must be either set to :monte_carlo or :historical"
        end

        if simulation.replications <= 0 &&
           (simulation.sim_type == :monte_carlo || simulation.randomize_years)
            errors *= "\nYou must specify a positive number of 'replications' for Monte Carlo or randomized historical simulations"
        end

        if simulation.number_of_cycles > 1 && simulation.reset_starting_levels == false
            errors *= "\nYou must not specify the 'number_of_cycles' if you do not reset the starting state values"
        end

        if simulation.number_of_cycles <= 0
            errors *= "\n'number_of_cycles' must be >= 1"
        end

        if simulation.reset_starting_levels ∉ [:default, false, true]
            errors *= "\n'reset_starting_levels' must be a boolean"
        end

        if simulation.initial_stage < 1 || simulation.initial_stage > WEEKSPERYEAR #change to number_of_wks
            errors *= "\n'initial_stage' must be between 1 and " * string(WEEKSPERYEAR)
        end

        if simulation.sim_type == :historical && simulation.sim_years == nothing
            errors *= "\n'sim_years' must be set for historical simulations"
        end

        if simulation.sim_type == :historical &&
           simulation.sim_years !== nothing &&
           length(simulation.sim_years) != simulation.replications
            errors *= "\n'replications' in a historical simulation must match the number of 'sim_years'"
        end
    end

    if rundata != nothing && simulation != nothing
        if !rundata.steady_state
            if simulation.number_of_cycles > 1
                errors *= "\nYou cannot simulate a finite-horizon model for more than 1 cycle."
            elseif simulation.reset_starting_levels == false
                errors *= "\nYou cannot simulate a finite-horizon model without reseting the starting levels."
            end
        end
    end

    if errors != ""
        @error("Invalid settings found:" * errors)
        error("Invalid settings found. See REPL for details.")
    end
end

"""
    get_file_directory(x::String, rundata::RunData; verbose::Bool = true)

This function looks in various locations for input file `x`. If the rundata references
a `scenario_dir` then the first directory to search is <data_dir>/data_files/<scenario_dir>.

### Required Arguments

`x` is the name of the input file.
`rundata` is a `JADE` `RunData` object.

### Keyword Arguments

`verbose` if set to `true` the function will print info to the REPL.
"""
function get_file_directory(x::String, rundata::RunData; verbose::Bool = true)
    input_directory = joinpath(@__JADE_DIR__, "Input", rundata.data_dir)

    if rundata.scenario_dir != "" &&
       !isdir(joinpath(input_directory, "data_files", rundata.scenario_dir))
        error(
            "'Input/$(rundata.data_dir)/data_files/$(rundata.scenario_dir)' scenario directory not found.",
        )
    elseif !isdir(joinpath(input_directory, "data_files"))
        if rundata.scenario_dir != ""
            error(
                "Cannot load scenario data files. 'Input/$(rundata.data_dir)/data_files' directory not found.",
            )
        end
        return joinpath(input_directory, x)
    elseif rundata.scenario_dir == ""
        return joinpath(input_directory, "data_files", x)
    elseif isfile(joinpath(input_directory, "data_files", rundata.scenario_dir, x))
        if verbose
            @info(
                "Loading scenario input file $x from $(joinpath(input_directory, "data_files", rundata.scenario_dir))"
            )
        end
        return joinpath(input_directory, "data_files", rundata.scenario_dir, x)
    else
        return joinpath(input_directory, "data_files", x)
    end
end
