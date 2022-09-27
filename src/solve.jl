#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
    optimize_policy!(
        JADEmodel::JADEModel,
        solveoptions::JADESolveOptions;
        async::Bool = false,
        print_level::Int = 1,
    )

This function trains the `JADEmodel`.

The output from this simulation will automatically be written into the
Output/`data_dir`/`policy_dir`/`sim_dir` subdirectory.

### Required Arguments
`JADEmodel` is the model detailing the year that we wish to create a policy for,
and all the corresponding data.

`solveoptions` contains all the SDDP.jl settings, as well as settings and
directories for loading cuts.

### Keyword Arguments
`async` is a boolean that sets whether SDDP.jl will run in parallel. Currently
recommended to be set to false.

`print_level` is an SDDP.jl setting.
"""
function optimize_policy!(
    JADEmodel::JADEModel,
    solveoptions::JADESolveOptions;
    async::Bool = false,
    print_level::Int = 1,
)
    d = JADEmodel.d
    sddpm = JADEmodel.sddpm
    previous_rundata = nothing

    check_settings_compatibility(rundata = d.rundata, solveoptions = solveoptions)

    if solveoptions.savedcuts != ""
        cuts_path = joinpath(@JADE_DIR, "Output", solveoptions.savedcuts, "cuts.json")
        if isfile(cuts_path) && solveoptions.warmstart_cuts
            previous_rundata = load_model_parameters(solveoptions.savedcuts)
        end
    else
        cuts_path = joinpath(
            @JADE_DIR,
            "Output",
            d.rundata.data_dir,
            d.rundata.policy_dir,
            "cuts.json",
        )
        if isfile(cuts_path) && solveoptions.warmstart_cuts
            previous_rundata =
                load_model_parameters(d.rundata.data_dir, d.rundata.policy_dir)
        end
    end

    @info(
        "Optimizing policy '" *
        d.rundata.policy_dir *
        "' for model data in " *
        joinpath("Input", d.rundata.data_dir)
    )

    if solveoptions.eoh_cutfile != "" && !solveoptions.warmstart_cuts
        final_week = mod(d.rundata.start_wk + d.rundata.number_of_wks - 2, 52) + 1

        cuts_path = joinpath(
            @JADE_DIR,
            "Input",
            d.rundata.data_dir,
            "EOH",
            "$(solveoptions.eoh_cutfile).eoh",
        )

        if !isfile(cuts_path)
            error("$(d.rundata.data_dir)/EOH/$(solveoptions.eoh_cutfile).eoh not found")
        end

        @info("Loading week " * string(final_week) * " cuts from " * cuts_path)

        previous_rundata = load_model_parameters(
            joinpath(
                @JADE_DIR,
                "Input",
                d.rundata.data_dir,
                "EOH",
                "$(solveoptions.eoh_cutfile).rdt",
            ),
        )

        check_rundata(d.rundata, previous_rundata, :eoh)

        JADE.read_finalcuts_from_file(sddpm, cuts_path, final_week, d.rundata.number_of_wks)
    elseif isfile(cuts_path)
        if solveoptions.warmstart_cuts
            if length(sddpm.nodes[1].bellman_function.global_theta.cuts) == 0 ||
               solveoptions.loadcuts_if_nonempty
                @info("Loading cuts from " * cuts_path)
                check_rundata(d.rundata, previous_rundata, :full)
                SDDP.read_cuts_from_file(sddpm, cuts_path)
                if has_upper_bound(
                    sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta,
                )
                    delete_upper_bound(
                        sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta,
                    )
                end
            else
                error(
                    "By default 'warmstart_cuts' only loads cuts into a model without previous cuts.\n Set 'loadcuts_if_nonempty = true' to override.",
                )
            end
        else
            if length(sddpm.nodes[1].bellman_function.global_theta.cuts) != 0
                @info("Existing cuts detected in model; these cuts will be retained")
            end
            cuts_path = joinpath(
                @JADE_DIR,
                "Output",
                d.rundata.data_dir,
                d.rundata.policy_dir,
                "cuts.json",
            )
            @info("Clearing previous cuts file")
            rm(cuts_path)
        end
    else
        if solveoptions.warmstart_cuts
            @info("No previous cuts file found")
        end
        if length(sddpm.nodes[1].bellman_function.global_theta.cuts) != 0
            @info("Existing cuts detected in model; these cuts will be retained")
        end
    end

    # Set seed
    Random.seed!(solveoptions.seed)

    # Policy generation
    if solveoptions.iterations > 0
        inflow_file = get_file_directory("inflows.csv", d.rundata)

        # Get our inflows, adjusted using DIA
        adjusted_inflows, firstweekinflows = adjustinflows(inflow_file, d.rundata)
        # Get inflows into the array form we will use in the JADE model
        inflow_mat = getinflows(adjusted_inflows, firstweekinflows, d.rundata)

        hist_inflow_mat = nothing
        sequences = nothing
        if solveoptions.fractionMC != 1.0
            sequences = Vector{Vector{Int}}
            seq_path = joinpath(
                @JADE_DIR,
                "Input",
                d.rundata.data_dir,
                solveoptions.custom_inflow_file,
            )
            if solveoptions.fractionMC < 1.0 && isfile(seq_path)
                sequences = getSequences(seq_path)
                @info(
                    "Loaded " *
                    string(length(sequences)) *
                    " custom inflow sequences for forward passes " *
                    "from " *
                    solveoptions.custom_inflow_file
                )
            else
                error(
                    "Custom inflow sequence file: " *
                    solveoptions.custom_inflow_file *
                    " not found.",
                )
            end
            if sequences != nothing
                if length(sequences[1]) != d.rundata.number_of_wks
                    error(
                        "Length of custom inflow sequences must match the number of weeks being modelled.",
                    )
                end
                min_year = minimum(minimum.(sequences))
                max_year = maximum(maximum.(sequences))
                hist_inflow_mat = getinflows_for_historical(
                    inflow_file,
                    d.rundata,
                    collect(min_year:max_year),
                )
            end
        else
            @info("No custom inflow sequences loaded for forward passes")
        end

        sample_paths = Vector{Tuple{Int,Dict{Symbol,Float64}}}[]
        count = 0
        extra = (d.rundata.steady_state && !solveoptions.reset_starting_levels) ? 1 : 0

        if d.rundata.steady_state && solveoptions.reset_starting_levels
            backup = sddpm.nodes[d.rundata.number_of_wks].children
            sddpm.nodes[d.rundata.number_of_wks].children = SDDP.Noise{Int64}[]
        end

        @info("Generating a policy...")
        for i in 1:solveoptions.iterations
            sample_path = Tuple{Int,Dict{Symbol,Float64}}[]
            if sequences == nothing || Random.rand() < solveoptions.fractionMC
                method = :montecarlo
                rand_years = rand(1:d.rundata.nscenarios, d.rundata.number_of_wks + extra)
            else
                method = :custom
                count = count % (length(sequences)) + 1
            end
            for t in 1:d.rundata.number_of_wks+extra
                s_inflows = Dict{Symbol,Float64}()
                if method == :montecarlo
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        s_inflows[c] =
                            inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][rand_years[t]]
                    end
                elseif method == :custom
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        if t <= d.rundata.number_of_wks
                            s_inflows[c] =
                                hist_inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][sequences[count][t]-min_year+1]
                        else
                            s_inflows[c] =
                                hist_inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][1]
                        end
                    end
                end
                push!(sample_path, ((t - 1) % d.rundata.number_of_wks + 1, s_inflows))
            end
            push!(sample_paths, sample_path)
        end

        parallel_scheme = nothing
        if async
            if d.parallel_optimizer == nothing
                parallel_scheme = SDDP.Asynchronous()
            else
                parallel_scheme = SDDP.Asynchronous() do m
                    optimizer = d.parallel_optimizer()
                    for node in values(m.nodes)
                        set_optimizer(node.subproblem, optimizer)
                    end
                end
            end
        else
            parallel_scheme = SDDP.Serial()
        end

        if d.rundata.steady_state && !solveoptions.reset_starting_levels
            solveresults = SDDP.train(
                sddpm,
                iteration_limit = solveoptions.iterations,
                cut_deletion_minimum = solveoptions.cutselection,
                sampling_scheme = WrapHistorical(sample_paths),
                cycle_discretization_delta = 10.0,
                dashboard = false,
                risk_measure = solveoptions.riskmeasure,
                parallel_scheme = parallel_scheme,
                print_level = print_level,
                forward_pass = JADEForwardPass(),
            )
        else
            solveresults = SDDP.train(
                sddpm,
                iteration_limit = solveoptions.iterations,
                cut_deletion_minimum = solveoptions.cutselection,
                sampling_scheme = SDDP.Historical(sample_paths),
                dashboard = false,
                risk_measure = solveoptions.riskmeasure,
                parallel_scheme = parallel_scheme,
                print_level = print_level,
                forward_pass = JADEForwardPass(),
            )
        end
        # Save cuts to a file
        write_training_results(sddpm, d, solveoptions)
    end

    if d.rundata.steady_state && solveoptions.reset_starting_levels
        sddpm.nodes[d.rundata.number_of_wks].children = backup
    end

    return nothing
end
