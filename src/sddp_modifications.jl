#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
	function read_finalcuts_from_file(
		model::SDDP.PolicyGraph{T},
		filename::String,
		cutsfromstage::Int,
		laststage::Int;
		node_name_parser::Function = SDDP._node_name_parser,
	) where {T}

### Required Arguments

`model` The SDDP `PolicyGraph` that cuts are being loaded into.

`filename` the full path to the cuts file being read.

`cutsfromstage` the index of the stage (not necessarily the week of the year) that cuts should be extracted from.

`laststage` the index of the stage of the model that the cuts are being loaded into.

`node_name_parser` SDDP function to convert the `Int` for laststage into the name of the SDDP node.

"""
function read_finalcuts_from_file(
    model::SDDP.PolicyGraph{T},
    filename::String,
    cutsfromstage::Int,
    laststage::Int;
    node_name_parser::Function = SDDP._node_name_parser,
) where {T}
    node_name = node_name_parser(T, string(laststage))::T
    node = model[node_name]
    bf = node.bellman_function
    if JuMP.has_upper_bound(bf.global_theta.theta)
        JuMP.delete_upper_bound(bf.global_theta.theta)
    end
    cuts = JSON.parsefile(filename, use_mmap = false)
    for node_cuts in cuts
        if node_cuts["node"] != string(cutsfromstage)
            continue
        end
        # Loop through and add the single-cuts.
        for json_cut in node_cuts["single_cuts"]
            has_state = haskey(json_cut, "state")
            state = if has_state
                Dict(Symbol(k) => v for (k, v) in json_cut["state"])
            else
                Dict(Symbol(k) => 0.0 for k in keys(json_cut["coefficients"]))
            end
            SDDP._add_cut(
                bf.global_theta,
                json_cut["intercept"],
                Dict(Symbol(k) => v for (k, v) in json_cut["coefficients"]),
                state,
                nothing,
                nothing;
                cut_selection = has_state,
            )
        end
        # Loop through and add the multi-cuts. There are two parts:
        #  (i) the cuts w.r.t. the state variable x
        # (ii) the cuts that define the risk set
        # There is one additional complication: if these cuts are being read
        # into a new model, the local theta variables may not exist yet.
        if length(node_cuts["risk_set_cuts"]) > 0
            SDDP._add_locals_if_necessary(
                node,
                bf,
                length(first(node_cuts["risk_set_cuts"])),
            )
        end
        for json_cut in node_cuts["multi_cuts"]
            has_state = haskey(json_cut, "state")
            state = if has_state
                Dict(Symbol(k) => v for (k, v) in json_cut["state"])
            else
                Dict(Symbol(k) => 0.0 for k in keys(json_cut["coefficients"]))
            end
            SDDP._add_cut(
                bf.local_thetas[json_cut["realization"]],
                json_cut["intercept"],
                Dict(Symbol(k) => v for (k, v) in json_cut["coefficients"]),
                state,
                nothing,
                nothing;
                cut_selection = has_state,
            )
        end
        # Here is part (ii): adding the constraints that define the risk-set
        # representation of the risk measure.
        for json_cut in node_cuts["risk_set_cuts"]
            expr = JuMP.@expression(
                node.subproblem,
                bf.global_theta.theta -
                sum(p * V.theta for (p, V) in zip(json_cut, bf.local_thetas))
            )
            if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
                JuMP.@constraint(node.subproblem, expr >= 0)
            else
                JuMP.@constraint(node.subproblem, expr <= 0)
            end
        end
    end
    return
end

"""
	function write_cuts_to_file(
		model::SDDP.PolicyGraph{T},
		filename::String,
		start_wk::Int,
	) where {T}

This function writes a cuts json file that is modified so that stage 1 is always associated with week 1.
This can only be used with a 52-week model. (It should only be used with steady-state JADE models.)

### Required Arguments

`model` The SDDP `PolicyGraph` that cuts are being written from.

`filename` the full path to the cuts file being written.

`start_wk` the week of the year that corresponds to stage 1 in the SDDP model.
"""
function write_cuts_to_file(
    model::SDDP.PolicyGraph{T},
    filename::String,
    start_wk::Int,
) where {T}
    if length(model.nodes) != 52
        error("This version of write_cuts_to_file must be used with a 52-week model.")
    end
    cuts = Dict{String,Any}[]
    for (node_name, node) in model.nodes
        if node.objective_state !== nothing || node.belief_state !== nothing
            error(
                "Unable to write cuts to file because model contains " *
                "objective states or belief states.",
            )
        end
        node_cuts = Dict(
            "node" => string(mod(node_name + start_wk - 2, 52) + 1),
            "single_cuts" => Dict{String,Any}[],
            "multi_cuts" => Dict{String,Any}[],
            "risk_set_cuts" => Vector{Float64}[],
        )
        oracle = node.bellman_function.global_theta
        for (cut, state) in zip(oracle.cuts, oracle.sampled_states)
            intercept = cut.intercept
            for (key, π) in cut.coefficients
                intercept += π * state.state[key]
            end
            push!(
                node_cuts["single_cuts"],
                Dict(
                    "intercept" => intercept,
                    "coefficients" => copy(cut.coefficients),
                    "state" => copy(state.state),
                ),
            )
        end
        for (i, theta) in enumerate(node.bellman_function.local_thetas)
            for (cut, state) in zip(theta.cuts, theta.sampled_states)
                intercept = cut.intercept
                for (key, π) in cut.coefficients
                    intercept += π * state.state[key]
                end
                push!(
                    node_cuts["multi_cuts"],
                    Dict(
                        "realization" => i,
                        "intercept" => intercept,
                        "coefficients" => copy(cut.coefficients),
                        "state" => copy(state.state),
                    ),
                )
            end
        end
        for p in node.bellman_function.risk_set_cuts
            push!(node_cuts["risk_set_cuts"], p)
        end
        push!(cuts, node_cuts)
    end
    open(filename, "w") do io
        return write(io, JSON.json(cuts))
    end
    return
end
