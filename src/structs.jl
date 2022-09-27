#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#-----------------------------------------------------
# Solving options
#-----------------------------------------------------

"""
An object containing SDDP settings used when training a JADE model.

### Fields

`iterations::Int` sddp iterations to be performed.

`riskmeasure::SDDP.AbstractRiskMeasure` sddp risk measure to use (set up to be CVaR).

`warmstart_cuts::Bool` if true will load existing cuts, either from local directory or 'savedcuts'.

`reset_starting_levels::Union{Symbol,Bool}` if true will reset starting storage levels to initial values each forward pass.

`seed::Int` random seed.

`savedcuts::String` alternate folder from which to retrieve cuts.

`cutselection::Int` sddp `cut_deletion_minimum` parameter.

`fractionMC::Float64` fraction of forward passes using monte-carlo.

`custom_inflow_file::String` file specifying sequences of inflows to be used in forward pass (fractionMC must be < 1).

`loadcuts_if_nonempty::Bool` set this to true if the cut file should be loaded into a model with existing cuts.

`eoh_cutfile::String` EOH cut file to load as the cuts for the final stage of a finite-horizon model.

`write_eohcuts::Bool` if true a steady-state model will write a set of files to Input/<data_dor>/EOH.
"""
mutable struct JADESolveOptions
    iterations::Int
    riskmeasure::SDDP.AbstractRiskMeasure
    warmstart_cuts::Bool
    reset_starting_levels::Union{Symbol,Bool}
    seed::Int
    savedcuts::String
    cutselection::Int
    fractionMC::Float64
    custom_inflow_file::String
    loadcuts_if_nonempty::Bool
    eoh_cutfile::String
    write_eohcuts::Bool
end

struct DecisionRule
    slope::Real
    intercept::Real
    boundtype::Symbol
    station::Symbol
    reservoir::Symbol
    flowtype::Symbol
    weeks::Vector{Int}
end

"""
    function DecisionRule(
        slope::Real,
        intercept::Real,
        boundtype::Symbol,
        station::Symbol,
        reservoir::Symbol,
        flowtype::Symbol,
        weeks::T where {T<:Union{Tuple{Int,Int},Int,UnitRange{Int},Vector}} = 1:52,
    )

This function constructs a DecisionRule object. See the JADE user documentation for
further details. DecisionRules can be applied to spill and generation at hydro stations,
these are both measured in MWh.

### Required Arguments
`slope` is the slope, w.r.t. the current storage in `reservoir`, of the (linear) bound applied
to generation and / or spill flow corresponding to `station`. Units are MWh / Mm³.

`intercept` is the bound applied to generation and / or spill flow corresponding to `station`
when `reservoir` is at 0 storage.  Units are (MWh).

`boundtype` is set to `:upper`, `:lower`, or `:equality` depending on how the flow should be
restricted.

`station` is the hydro station that the decision rule restricts.

`reservoir` is the hydro reservoir whose  affects the bound applied.

`flowtype` is set to `:generation`, `:spill`, or `:combined`.

`weeks` defines the weeks in the year that the decision rule is in effect. E.g. [1:10,43:52].
"""
function DecisionRule(
    slope::Real,
    intercept::Real,
    boundtype::Symbol,
    station::Symbol,
    reservoir::Symbol,
    flowtype::Symbol,
    weeks::T where {T<:Union{Tuple{Int,Int},Int,UnitRange{Int},Vector}} = 1:52,
)
    if boundtype ∉ [:lower, :upper, :equality]
        error("Invalid bound type; it must be :lower, :upper or :equality.")
    end

    if flowtype ∉ [:generation, :spill, :combined]
        error("Invalid flow type; it must be :generation, :spill or :combined.")
    end

    return DecisionRule(
        slope,
        intercept,
        boundtype,
        station,
        reservoir,
        flowtype,
        convert_to_int_array(weeks),
    )
end

#-----------------------------------------------------
# Run-related information
#-----------------------------------------------------
"""
An object containing the JADE model settings, required to define the stage problems in SDDP.jl.

### Fields

`data_dir::String` directory for writing output.

`start_yr::Int` first year in the problem.

`start_wk::Int` first week in the problem.

`number_of_wks::Int` the number of weeks that we solve for.

`sample_years::Vector{Int}` array of years that inflows are sampled from.

`dialength::Int` inflow correlation length.

`penalty_lb::Float64` penalty for flows below lower bound.

`penalty_ub::Float64` penalty for flows above upper bound.

`policy_dir::String` output directory for that the policy cuts and settings are stored in.

`discount::Float64` steady-state discount factor (0 if finite-horizon).

`losses::Symbol` type of losses (:none,:piecewise,:resistive).

`nscenarios::Int` number of possible inflow realisations per week.

`weekly_discounting::Bool` set to true if (equivalent) discount should applied weekly.

`scale_reservoirs::Float64` scale factor for reservoirs.

`scale_objective::Float64` scale factor for objective.

`use_terminal_mwvs::Bool` boolean specifying whether we use the given terminal MWVs.

`first_week_known::Bool` set to true if first-week inflow is set to historical value from training year.

`steady_state::Bool` set to true if the model is steady-state.

`decision_rules::Vector{DecisionRule}` vector of decision rules that can restrict how reservoirs are operated.

`scenario_dir::String` name of directory within Input/<data_dir>/data_files/ containing data files that change within the scenario.
"""
mutable struct RunData
    data_dir::String
    start_yr::Int
    start_wk::Int
    number_of_wks::Int
    sample_years::Vector{Int}
    dialength::Int
    penalty_lb::Float64
    penalty_ub::Float64
    policy_dir::String
    discount::Float64
    losses::Symbol
    nscenarios::Int
    weekly_discounting::Bool
    scale_reservoirs::Float64
    scale_objective::Float64
    use_terminal_mwvs::Bool
    first_week_known::Bool
    steady_state::Bool
    decision_rules::Vector{DecisionRule}
    scenario_dir::String
end

mutable struct Tranche
    q::Float64
    p::Float64
end

mutable struct Sets
    BLOCKS::Vector{Symbol}                  # load blocks
    NODES::Vector{Symbol}                   # power system nodes
    SECTORS::Vector{Symbol}                 # load sectors
    THERMALS::Vector{Symbol}                # thermal power stations
    RESERVOIRS::Vector{Symbol}              # reservoirs (can store water)
    JUNCTIONS::Vector{Symbol}               # hydro junctions
    CATCHMENTS::Vector{Symbol}              # all junctions and reservoirs
    HYDROS::Vector{Symbol}                  # hydro power stations
    CATCHMENTS_WITH_INFLOW::Vector{Symbol}  # locations with inflow data
    JUNCTIONS_WITHOUT_INFLOW::Vector{Symbol}
    NATURAL_ARCS::Vector{NTuple{2,Symbol}} # arcs independent of hydro stations
    STATION_ARCS::Vector{NTuple{2,Symbol}} # origin and destination of water for a hydro station
    TRANS_ARCS::Vector{NTuple{2,Symbol}}   # power transmission arcs
end

function Sets()
    return Sets(
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        Symbol[],
        NTuple{2,Symbol}[],
        NTuple{2,Symbol}[],
        NTuple{2,Symbol}[],
    )
end

"""
An object containing JADE simulation settings.

### Fields

`sim_dir::String` is the folder that the simulation output should be written to (`Output/<data_dir>/<policy_dir>/<sim_dir>`).

`sim_type::Symbol` is the type of simulation, either `:monte_carlo` or `:historical`.

`replications::Int` is the number of replications of the simulation.

`sim_years::Union{Nothing,Vector{Int}}` is either `nothing` or a vector of years that should be sampled from in the simulation.

`randomize_years::Bool` if `true` the simulation will randomly choose years from `sim_years` (with replacement); otherwise the simulation sequentially samples the years from `sim_years` (in this case replications is automatically set).

`reset_starting_levels::Union{Symbol,Bool}` if `true` the simulation will reset the initial storage levels after each iteration, even for steady-state models.

`number_of_cycles::Int` for steady-state models, if this is greater than 1, each replication of the simulation will run `number_of_cycles` times; `reset_starting_levels` must be set to `true`.

`initial_stage::Int` is the stage (not the week) of the simulation that the simulation will start from.

`initial_state::Union{Nothing,Dict{String,Float64}}` is the initial state variables in the `initial stage`

`random_seed::Int` is the random seed for the simulation.
"""
mutable struct JADESimulation
    sim_dir::String
    sim_type::Symbol
    replications::Int
    sim_years::Union{Nothing,Vector{Int}}
    randomize_years::Bool
    reset_starting_levels::Union{Symbol,Bool}
    number_of_cycles::Int
    initial_stage::Int
    initial_state::Union{Nothing,Dict{String,Float64}}
    random_seed::Int
end
