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
