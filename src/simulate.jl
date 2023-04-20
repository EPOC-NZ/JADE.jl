#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
	simulate(JADEmodel::JADEModel, parameters::JADESimulation)

This function carries out a simulation of a JADE model; all the specifications
for this simulation are defined within the `parameters` argument.

The output from this simulation will automatically be written into the
Output/`data_dir`/`policy_dir`/`sim_dir` subdirectory.

### Required Arguments
`JADEmodel` is the model detailing the scenario that we wish to simulate, and
all the corresponding data.

`parameters` contains all the simulation information, including the number of
replications, the type of simulation, the hydrological years to sample from, etc.
"""
function simulate(JADEmodel::JADEModel, parameters::JADESimulation)
    d = JADEmodel.d
    sddpm = JADEmodel.sddpm

    check_settings_compatibility(rundata = d.rundata, simulation = parameters)

    @info(
        "Simulating policy '" *
        d.rundata.policy_dir *
        "' for model data in " *
        joinpath("Input", d.rundata.data_dir)
    )

    cuts_path = joinpath(
        @__JADE_DIR__,
        "Output",
        d.rundata.data_dir,
        d.rundata.policy_dir,
        "cuts.json",
    )

    if length(sddpm.nodes[1].bellman_function.global_theta.cuts) == 0
        if isfile(cuts_path)
            @info(
                "Loading existing cuts file: " *
                joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir, "cuts.json")
            )
            previous_rundata =
                load_model_parameters(d.rundata.data_dir, d.rundata.policy_dir)
            check_rundata(d.rundata, previous_rundata, :partial)
            SDDP.read_cuts_from_file(sddpm, cuts_path)
            if JuMP.has_upper_bound(
                sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta,
            )
                JuMP.delete_upper_bound(
                    sddpm.nodes[d.rundata.number_of_wks].bellman_function.global_theta.theta,
                )
            end
        else
            @info(
                "No cuts file found in " *
                joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir)
            )
        end
    end

    get_primal = [
        :thermal_use,
        :hydro_disp,
        :naturalflows,
        :releases,
        :spills,
        :reslevel,
        :inflow,
        :flowover,
        :flowunder,
        :flowpenalties,
        :lostload,
        :terminalcost,
        :transflow,
        :lostloadcosts,
        :contingent_storage_cost,
        :carbon_emissions,
    ]

    get_dual = Dict{Symbol,Function}(
        :prices => (sp) -> d.rundata.scale_objective * JuMP.dual.(sp[:defineShedding]),
        :mwv =>
            (sp) ->
                -d.rundata.scale_objective * JuMP.dual.(sp[:rbalance]) / 1E3 /
                d.rundata.scale_reservoirs,
    )

    initial_state = Dict(String(k) => v for (k, v) in JADEmodel.sddpm.initial_root_state)
    if parameters.initial_state == nothing
        @info("Using default initial reservoir levels")
    else
        if length(initial_state) != length(parameters.initial_state)
            @info("An initial value must be specified for each of the following states:")
            for key in keys(initial_state)
                @info(key[10:end-1])
            end
            error("initial_state dictionary has incorrect number of states")
        end
        for key in keys(parameters.initial_state)
            if key ∉ keys(initial_state)
                @info(
                    "An initial value must be specified for each of the following states:"
                )
                for key in keys(initial_state)
                    @info(key[10:end-1])
                end
                error("Invalid state: " * key)
            end
        end
        initial_state = parameters.initial_state
        @info("Using simulation-specific initial reservoir levels")
    end

    if parameters.sim_type == :monte_carlo
        @info(
            "Simulating Monte Carlo inflows with " *
            string(parameters.replications) *
            " replications"
        )
    elseif parameters.sim_type == :historical
        if parameters.randomize_years == false
            @info("Simulating historical inflows in sequence provided")
        else
            @info(
                "Simulating annual historical inflows randomly with " *
                string(parameters.replications) *
                " replications"
            )
        end
    end

    Random.seed!(parameters.random_seed)
    wks = d.rundata.number_of_wks * parameters.number_of_cycles
    lastwk = d.rundata.number_of_wks * parameters.number_of_cycles
    if !d.rundata.steady_state
        wks += 1 - parameters.initial_stage
    else
        lastwk += parameters.initial_stage - 1
    end

    if parameters.sim_type == :monte_carlo
        if !d.rundata.steady_state || parameters.reset_starting_levels == true
            results = SDDP.simulate(
                sddpm,
                parameters.replications,
                get_primal,
                custom_recorders = get_dual,
                sampling_scheme = SDDP.InSampleMonteCarlo(
                    max_depth = wks,
                    initial_node = parameters.initial_stage,
                    terminate_on_cycle = false,
                    terminate_on_dummy_leaf = false,
                ),
                incoming_state = initial_state,
            )

            for i in 1:parameters.replications
                for τ in parameters.initial_stage:lastwk
                    t = τ - parameters.initial_stage + 1
                    for key in keys(results[i][t][:reslevel])
                        results[i][t][:reslevel][key[1]] = SDDP.State(
                            results[i][t][:reslevel][key[1]].in *
                            d.rundata.scale_reservoirs,
                            results[i][t][:reslevel][key[1]].out *
                            d.rundata.scale_reservoirs,
                        )
                    end
                    if results[i][t][:noise_term][:scenario] == 0
                        results[i][t][:inflow_year] = d.rundata.start_yr
                    else
                        results[i][t][:inflow_year] = d.rundata.sample_years[round(
                            Int,
                            results[i][t][:noise_term][:scenario],
                        )]
                    end
                end
            end
        else
            sequence = SDDP.simulate(
                sddpm,
                1,
                get_primal,
                custom_recorders = get_dual,
                sampling_scheme = SDDP.InSampleMonteCarlo(
                    max_depth = wks * parameters.replications,
                    initial_node = parameters.initial_stage,
                    terminate_on_cycle = false,
                    terminate_on_dummy_leaf = false,
                ),
                incoming_state = initial_state,
            )

            results = Vector{Dict{Symbol,Any}}[]
            for i in 1:parameters.replications
                push!(results, Dict{Symbol,Any}[])
                for τ in parameters.initial_stage:lastwk
                    t = τ - parameters.initial_stage + 1
                    temp = Dict{Symbol,Any}()
                    push!(results[i], temp)
                    for sym in
                        vcat(get_primal, [:stage_objective, :bellman_term, :prices, :mwv])
                        if sym == :reslevel
                            results[i][t][sym] = Dict{Symbol,SDDP.State}()
                            for key in keys(sequence[1][(i-1)*wks+t][sym])
                                results[i][t][sym][key[1]] = SDDP.State(
                                    sequence[1][(i-1)*wks+t][sym][key].in *
                                    d.rundata.scale_reservoirs,
                                    sequence[1][(i-1)*wks+t][sym][key].out *
                                    d.rundata.scale_reservoirs,
                                )
                            end
                        else
                            results[i][t][sym] = sequence[1][(i-1)*wks+t][sym]
                        end
                    end
                    if sequence[1][(i-1)*wks+t][:noise_term][:scenario] == 0
                        results[i][t][:inflow_year] = d.rundata.start_yr
                    else
                        results[i][t][:inflow_year] = d.rundata.sample_years[round(
                            Int,
                            sequence[1][(i-1)*wks+t][:noise_term][:scenario],
                        )]
                    end
                end
            end
        end
    elseif parameters.sim_type == :historical
        sequence = nothing
        results = Vector{Dict{Symbol,Any}}[]
        if !d.rundata.steady_state || parameters.reset_starting_levels == true
            sample_paths = Vector{Tuple{Int,Dict{Symbol,Float64}}}[]
            push!(sample_paths, Tuple{Int,Dict{Symbol,Float64}}[])
            seq = []
            count = 0
            for year in parameters.sim_years
                years = collect(
                    year:(year-1+ceil(
                        Int,
                        (d.rundata.start_wk - 1 + lastwk) / WEEKSPERYEAR,
                    )),
                )

                inflow_mat = getinflows_for_historical(
                    get_file_directory("inflows.csv", d.rundata),
                    d.rundata,
                    years,
                )
                i = 1
                extrawks = d.rundata.steady_state ? parameters.initial_stage - 1 : 0
                for t in parameters.initial_stage:(d.rundata.number_of_wks+extrawks)
                    s_inflows = Dict{Symbol,Float64}()
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        s_inflows[c] =
                            inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][i]
                    end
                    if d.rundata.first_week_known && t == 1
                        s_inflows[:scenario] = d.rundata.start_yr
                    else
                        s_inflows[:scenario] = years[i]
                    end
                    if (t + d.rundata.start_wk - 2) % WEEKSPERYEAR == WEEKSPERYEAR - 1
                        i += 1
                    end
                    push!(sample_paths[end], ((t - 1) % WEEKSPERYEAR + 1, s_inflows))
                end
                count += 1
                if count == parameters.number_of_cycles
                    count = 0
                    push!(sample_paths, Tuple{Int,Dict{Symbol,Float64}}[])
                end
            end
            sims = SDDP.simulate(
                sddpm,
                parameters.replications,
                get_primal,
                sampling_scheme = SDDP.Historical(sample_paths),
                custom_recorders = get_dual,
                incoming_state = initial_state,
            )

            for sim in sims
                append!(seq, sim)
            end
            sequence = [seq]
        else
            sample_path = Tuple{Int,Dict{Symbol,Float64}}[]
            for year in parameters.sim_years
                years = collect(
                    year:(year-1+ceil(
                        Int,
                        (
                            d.rundata.start_wk + d.rundata.number_of_wks - 1 +
                            parameters.initial_stage - 1
                        ) / WEEKSPERYEAR,
                    )),
                )

                inflow_mat = getinflows_for_historical(
                    get_file_directory("inflows.csv", d.rundata),
                    d.rundata,
                    years,
                )
                i = 1

                for t in
                    parameters.initial_stage:(d.rundata.number_of_wks+parameters.initial_stage-1)
                    s_inflows = Dict{Symbol,Float64}()
                    for c in d.sets.CATCHMENTS_WITH_INFLOW
                        s_inflows[c] =
                            inflow_mat[(t+d.rundata.start_wk-2)%WEEKSPERYEAR+1][c][i]
                    end
                    s_inflows[:scenario] = years[i]
                    if (t + d.rundata.start_wk - 2) % WEEKSPERYEAR == WEEKSPERYEAR - 1
                        i += 1
                    end
                    push!(sample_path, ((t - 1) % WEEKSPERYEAR + 1, s_inflows))
                end
            end
            sequence = SDDP.simulate(
                sddpm,
                1,
                get_primal,
                sampling_scheme = SDDP.Historical(sample_path),
                custom_recorders = get_dual,
                incoming_state = initial_state,
            )
        end

        for i in 1:parameters.replications
            push!(results, Dict{Symbol,Any}[])
            for τ in parameters.initial_stage:lastwk
                t = τ - parameters.initial_stage + 1
                temp = Dict{Symbol,Any}()
                push!(results[i], temp)
                for sym in
                    vcat(get_primal, [:stage_objective, :bellman_term, :prices, :mwv])
                    if sym == :reslevel
                        results[i][t][sym] = Dict{Symbol,SDDP.State}()
                        for key in keys(sequence[1][(i-1)*wks+t][sym])
                            results[i][t][sym][key[1]] = SDDP.State(
                                sequence[1][(i-1)*wks+t][sym][key].in *
                                d.rundata.scale_reservoirs,
                                sequence[1][(i-1)*wks+t][sym][key].out *
                                d.rundata.scale_reservoirs,
                            )
                        end
                    else
                        results[i][t][sym] = sequence[1][(i-1)*wks+t][sym]
                    end
                end
                results[i][t][:inflow_year] =
                    round(Int, sequence[1][(i-1)*wks+t][:noise_term][:scenario])
            end
        end
    end

    for i in 1:parameters.replications
        for τ in parameters.initial_stage:lastwk
            t = τ - parameters.initial_stage + 1
            results[i][t][:stage_objective] *= d.rundata.scale_objective
            results[i][t][:bellman_term] *= d.rundata.scale_objective
            for r in d.sets.RESERVOIRS
                results[i][t][:mwv][r] /= d.reservoirs[r].sp
            end
            if t == 1
                results[i][t][:running_cost] = results[i][t][:stage_objective]
            else
                results[i][t][:running_cost] =
                    results[i][t-1][:running_cost] + results[i][t][:stage_objective]
            end

            results[i][t][:total_storage] = 0
            for r in keys(d.reservoirs)
                results[i][t][:total_storage] +=
                    results[i][t][:reslevel][r].out * d.reservoirs[r].sp * 1000
            end
        end
    end

    @info(
        "Saving output in " *
        joinpath("Output", d.rundata.data_dir, d.rundata.policy_dir, parameters.sim_dir)
    )
    write_sim_results(results, d, parameters)
    output_tidy_results(
        results,
        d,
        parameters,
        variables = [
            :reslevel,
            :thermal_use,
            :transflow,
            :prices,
            :lostload,
            :flowover,
            :flowunder,
            :flowpenalties,
            :contingent_storage_cost,
            :carbon_emissions,
            :spills,
            :total_storage,
            :inflow_year,
            :mwv,
        ],
    )

    @info("Done.")
    return results
end
