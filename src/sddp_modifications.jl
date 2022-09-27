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
    if has_upper_bound(bf.global_theta.theta)
        delete_upper_bound(bf.global_theta.theta)
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
            expr = @expression(
                node.subproblem,
                bf.global_theta.theta -
                sum(p * V.theta for (p, V) in zip(json_cut, bf.local_thetas))
            )
            if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
                @constraint(node.subproblem, expr >= 0)
            else
                @constraint(node.subproblem, expr <= 0)
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

mutable struct WrapHistorical{T,S} <: SDDP.AbstractSamplingScheme
    scenarios::Vector{SDDP.Noise{Vector{Tuple{T,S}}}}
    sequential::Bool
    counter::Int
end

function Base.show(io::IO, h::WrapHistorical)
    print(
        io,
        "A Historical sampler with $(length(h.scenarios)) scenarios sampled ",
        h.sequential ? "sequentially." : "probabilistically.",
    )
    return
end

"""
    WrapHistorical(
        scenarios::Vector{Vector{Tuple{T,S}}},
        probability::Vector{Float64},
    ) where {T,S}

A sampling scheme that samples a scenario from the vector of scenarios
`scenarios` according to `probability`.

### Example

```julia
WrapHistorical(
    [
        [(1, 0.5), (2, 1.0), (3, 0.5)],
        [(1, 0.5), (2, 0.0), (3, 1.0)],
        [(1, 1.0), (2, 0.0), (3, 0.0)]
    ],
    [0.2, 0.5, 0.3],
)
```
"""
function WrapHistorical(
    scenarios::Vector{Vector{Tuple{T,S}}},
    probability::Vector{Float64},
) where {T,S}
    if !(sum(probability) ≈ 1.0)
        error(
            "Probability of historical scenarios must sum to 1. Currently: " *
            "$(sum(probability)).",
        )
    end
    output = [SDDP.Noise(s, p) for (s, p) in zip(scenarios, probability)]
    return WrapHistorical(output, false, 0)
end

"""
    WrapHistorical(scenarios::Vector{Vector{Tuple{T,S}}}) where {T,S}

A deterministic sampling scheme that iterates through the vector of provided
`scenarios`. The initial states in each iteration (other than the first) are
set to the ending state of the previous forward pass.

Used in the forward pass of the training to enable custom inflow sequences in steady-state models (see [solve.jl](https://github.com/EPOC-NZ/JADE/blob/master/src/solve.jl)).

## Example

```julia
WrapHistorical([
    [(1, 0.5), (2, 1.0), (3, 0.5)],
    [(1, 0.5), (2, 0.0), (3, 1.0)],
    [(1, 1.0), (2, 0.0), (3, 0.0)],
])
```
"""
function WrapHistorical(scenarios::Vector{Vector{Tuple{T,S}}}) where {T,S}
    return WrapHistorical(SDDP.Noise.(scenarios, NaN), true, 0)
end

WrapHistorical(scenario::Vector{Tuple{T,S}}) where {T,S} = WrapHistorical([scenario])

"""
	SDDP.sample_scenario(
		graph::SDDP.PolicyGraph{T},
		sampling_scheme::WrapHistorical{T,NoiseTerm};
		kwargs...,
	) where {T,NoiseTerm}
	
A method to sample a scenario using historical data, that returns `true`
for `terminated_due_to_cycle`. This is used with a custom forward pass: 
`JADEForwardPass`.
"""
function SDDP.sample_scenario(
    graph::SDDP.PolicyGraph{T},
    sampling_scheme::WrapHistorical{T,NoiseTerm};
    # Ignore the other kwargs because the user is giving
    # us the full scenario.
    kwargs...,
) where {T,NoiseTerm}
    if sampling_scheme.sequential
        sampling_scheme.counter += 1
        if sampling_scheme.counter > length(sampling_scheme.scenarios)
            sampling_scheme.counter = 1
        end
        return sampling_scheme.scenarios[sampling_scheme.counter].term, true
    end
    return SDDP.sample_noise(sampling_scheme.scenarios), true
end

struct JADEForwardPass <: SDDP.AbstractForwardPass end

"""
	SDDP.forward_pass(
		model::SDDP.PolicyGraph{T},
		options::SDDP.Options,
		::JADEForwardPass,
	) where {T}

This custom forward pass method for JADE enables historical sequences of inflows to
wrap back, giving new starting states for the next iteration. In order to fit with SDDP.jl,
if we wish to enable this functionality, there is an additional node given in the historical
scenario path. This is used to determine `final_node`. We are not using Markov states, so this code
is a bit more complicated than is currently necessary, but should work if Markov states were 
to be introduced.
"""
function SDDP.forward_pass(
    model::SDDP.PolicyGraph{T},
    options::SDDP.Options,
    ::JADEForwardPass,
) where {T}
    # First up, sample a scenario. Note that if a cycle is detected, this will
    # return the cycle node as well.
    SDDP.TimerOutputs.@timeit SDDP.SDDP_TIMER "sample_scenario" begin
        scenario_path, terminated_due_to_cycle =
            SDDP.sample_scenario(model, options.sampling_scheme)
    end

    if terminated_due_to_cycle
        final_node = scenario_path[end]
        scenario_path = scenario_path[1:end-1]
    end

    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol,Float64}[]

    # Our initial incoming state.
    incoming_state_value = copy(options.initial_state)

    # A cumulator for the stage-objectives.
    cumulative_value = 0.0

    # Iterate down the scenario.
    for (depth, (node_index, noise)) in enumerate(scenario_path)
        node = model[node_index]

        # ===== Begin: starting state for infinite horizon =====
        starting_states = options.starting_states[node_index]
        if length(starting_states) > 0
            # There is at least one other possible starting state. If our
            # incoming state is more than δ away from the other states, add it
            # as a possible starting state.
            if SDDP.distance(starting_states, incoming_state_value) >
               options.cycle_discretization_delta
                push!(starting_states, incoming_state_value)
            end
            incoming_state_value = splice!(starting_states, rand(1:length(starting_states)))
        end
        # ===== End: starting state for infinite horizon =====
        # Solve the subproblem, note that `require_duals = false`.
        SDDP.TimerOutputs.@timeit SDDP.SDDP_TIMER "solve_subproblem" begin
            subproblem_results = SDDP.solve_subproblem(
                model,
                node,
                incoming_state_value,
                noise,
                scenario_path[1:depth],
                duality_handler = nothing,
            )
        end
        # Cumulate the stage_objective.
        cumulative_value += subproblem_results.stage_objective

        # Set the outgoing state value as the incoming state value for the next
        # node.
        incoming_state_value = copy(subproblem_results.state)
        # Add the outgoing state variable to the list of states we have sampled
        # on this forward pass.
        push!(sampled_states, incoming_state_value)
    end
    if terminated_due_to_cycle
        # Get the last node in the scenario.
        final_node_index = final_node[1]
        # We terminated due to a cycle. Here is the list of possible starting
        # states for that node:
        starting_states = options.starting_states[final_node_index]
        # We also need the incoming state variable to the final node, which is
        # the outgoing state value of the last node:
        incoming_state_value = sampled_states[end]
        # If this incoming state value is more than δ away from another state,
        # add it to the list.
        if SDDP.distance(starting_states, incoming_state_value) >
           options.cycle_discretization_delta
            push!(starting_states, incoming_state_value)
        end
    end
    # ===== End: drop off starting state if terminated due to cycle =====
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        #objective_states = objective_states,
        #belief_states = belief_states,
        objective_states = NTuple{1,Float64}[],
        belief_states = Tuple{Int,Dict{T,Float64}}[],
        cumulative_value = cumulative_value,
    )
end

struct InSampleMonteCarlo2 <: SDDP.AbstractSamplingScheme
    max_depth::Int
    initial_stage::Int
    terminate_on_cycle::Bool
    terminate_on_dummy_leaf::Bool
    rollout_limit::Function
end

"""
	function InSampleMonteCarlo2(;
		max_depth::Int = 0,
		initial_stage::Int = 1,
		terminate_on_cycle::Bool = false,
		terminate_on_dummy_leaf::Bool = true,
		rollout_limit::Function = i -> typemax(Int),
	)

This function creates an instance of the `InSampleMonteCarlo2` `sampling_scheme`.
This is used for Monte Carlo simulations that may start from an arbitrary stage
(`initial_stage`).

### Keyword Arguments

`max_depth` is the maximum number of stages to simulate.

`initial_stage` is the stage number of the initial stage.

`terminate_on_cycle` if `true` the simulation will stop if it returns to a previously visited node.

`terminate_on_dummy_leaf` if `true` the simulation will may stop if it the probabilities of the noises don't sum to 1.

`rollout_limit` is a function that allows the number of simulated stages to vary with the training iteration number.
"""
function InSampleMonteCarlo2(;
    max_depth::Int = 0,
    initial_stage::Int = 1,
    terminate_on_cycle::Bool = false,
    terminate_on_dummy_leaf::Bool = true,
    rollout_limit::Function = i -> typemax(Int),
)
    if !terminate_on_cycle && !terminate_on_dummy_leaf && max_depth == 0
        error(
            "terminate_on_cycle and terminate_on_dummy_leaf cannot both be " *
            "false when max_depth=0.",
        )
    end
    new_rollout = let i = 0
        () -> (i += 1; rollout_limit(i))
    end
    return InSampleMonteCarlo2(
        max_depth,
        initial_stage,
        terminate_on_cycle,
        terminate_on_dummy_leaf,
        new_rollout,
    )
end

"""
	function SDDP.sample_scenario(
		graph::SDDP.PolicyGraph{T},
		sampling_scheme::InSampleMonteCarlo2,
	) where {T}

This function is a version of the SDDP.jl `SDDP.sample_scenario` function which
uses the `InSampleMonteCarlo2` `sampling_scheme`. This enables a simulation to start
from an arbitrary stage (`initial_stage`).

### Required Arguments

`graph` is an SDDP.jl `PolicyGraph`.

`sampling_scheme` is an instance of `InSampleMonteCarlo2` sampling scheme.
"""
function SDDP.sample_scenario(
    graph::SDDP.PolicyGraph{T},
    sampling_scheme::InSampleMonteCarlo2,
) where {T}
    max_depth = min(sampling_scheme.max_depth, sampling_scheme.rollout_limit())
    # Storage for our scenario. Each tuple is (node_index, noise.term).
    scenario_path = Tuple{T,Any}[]
    # We only use visited_nodes if terminate_on_cycle=true. Just initialize
    # anyway.
    visited_nodes = Set{T}()
    node_index = sampling_scheme.initial_stage
    while true
        node = graph[node_index]
        noise_terms = node.noise_terms
        children = node.children
        noise = SDDP.sample_noise(noise_terms)
        push!(scenario_path, (node_index, noise))
        # Termination conditions:
        if length(children) == 0
            # 1. Our node has no children, i.e., we are at a leaf node.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_cycle && node_index in visited_nodes
            # 2. terminate_on_cycle = true and we have detected a cycle.
            return scenario_path, true
        elseif 0 < sampling_scheme.max_depth <= length(scenario_path)
            # 3. max_depth > 0 and we have explored max_depth number of nodes.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_dummy_leaf &&
               rand() < 1 - sum(child.probability for child in children)
            # 4. we sample a "dummy" leaf node in the next step due to the
            # probability of the child nodes summing to less than one.
            return scenario_path, false
        end
        # We only need to store a list of visited nodes if we want to terminate
        # due to the presence of a cycle.
        if sampling_scheme.terminate_on_cycle
            push!(visited_nodes, node_index)
        end
        # Sample a new node to transition to.
        node_index = SDDP.sample_noise(children)::T
    end
    # Throw an error because we should never end up here.
    return error("Internal SDDP error: something went wrong sampling a scenario.")
end
