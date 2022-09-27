#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
    computefixed(fixed_file::String, outage_file::String, demand, T)

# Description

This function computes the fixed generation (MWh) at every node.

Currently the fixed stations only support constant output capacity,
although outages can be set to enable time-varying output.
"""
function computefixed(fixed_file::String, demand, outage_data, T)
    # {station,node,block}=>quantity (MWh)
    all_demand_reduction = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
    tp_demand_reduction = Dict{Tuple{Symbol,Symbol,Symbol},Dict{TimePoint,Float64}}()
    blocks = Symbol[]
    parsefile(fixed_file, true) do items
        @assert length(items) >= 5
        if lowercase(items[1]) == "station"
            for it in items[5:end]
                push!(blocks, str2sym(it))
            end
            return
        end

        station = str2sym(items[1])
        node = str2sym(items[2])

        if items[3] != "all" || items[4] != "all"
            tp = TimePoint(parse(Int, items[3]), parse(Int, items[4]))
            for (i, blk) in enumerate(blocks)
                if haskey(all_demand_reduction, (station, node, blk))
                    error(
                        "Fixed generator " *
                        string(station) *
                        " defined for both single TimePeriods and all periods.",
                    )
                end
                if !haskey(tp_demand_reduction, (station, node, blk))
                    tp_demand_reduction[(station, node, blk)] = Dict{TimePoint,Float64}()
                end
                if !haskey(tp_demand_reduction[(station, node, blk)], tp)
                    tp_demand_reduction[(station, node, blk)][tp] = 0.0
                end
                tp_demand_reduction[(station, node, blk)][tp] += parse(Float64, items[4+i])
            end
        else
            for (i, blk) in enumerate(blocks)
                if haskey(tp_demand_reduction, (station, node, blk))
                    error(
                        "Fixed generator " *
                        string(station) *
                        " defined for both single TimePeriods and all periods.",
                    )
                end
                if !haskey(all_demand_reduction, (station, node, blk))
                    all_demand_reduction[(station, node, blk)] = 0.0
                end
                all_demand_reduction[(station, node, blk)] += parse(Float64, items[4+i])
            end
        end
    end
    fixed = deepcopy(demand)
    for (tp, f) in fixed
        for (node, blk) in keys(f)
            f[(node, blk)] = 0.0
        end
        for ((station, node, blk), val) in all_demand_reduction
            if haskey(f, (node, blk))
                # Note: we account for outages in the fixed stations
                if (station, blk) ∈ keys(outage_data[tp])
                    f[(node, blk)] += (val - outage_data[tp][(station, blk)]) * T[tp][blk]
                else
                    f[(node, blk)] += val * T[tp][blk]
                end
            end
        end
    end
    for ((station, node, blk), vals) in tp_demand_reduction
        for (tp, val) in vals
            f = fixed[tp]
            if haskey(f, (node, blk))
                # Note: we account for outages in the fixed stations
                if (station, blk) ∈ keys(outage_data[tp])
                    f[(node, blk)] += (val - outage_data[tp][(station, blk)]) * T[tp][blk]
                else
                    f[(node, blk)] += val * T[tp][blk]
                end
            end
        end
    end
    return fixed
end

function getfixedgeneration(
    fixed_file::String,
    years::Vector{Int},
    weeks::Vector{Int},
    loadblocks::Vector{Symbol},
    demand::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    outage_data::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    T::TimeSeries{Dict{Symbol,Float64}},
)
    tp_demand_reduction = Dict{Tuple{Symbol,Symbol,Symbol},Dict{TimePoint,Float64}}()

    parsefile(fixed_file, true) do items
        @assert length(items) == 6
        if lowercase(items[1]) == "station"
            return
        end

        station = str2sym(items[1])
        node = str2sym(items[2])

        years_ = process_set_items(items[3], years)
        weeks_ = process_set_items(items[4], weeks)
        loadblocks_ = process_set_items(items[5], loadblocks)
        generation = parse(Float64, items[6])

        for blk in loadblocks_
            if !haskey(tp_demand_reduction, (station, node, blk))
                tp_demand_reduction[(station, node, blk)] = Dict{TimePoint,Float64}()
            end
            for y in years_
                for w in weeks_
                    tp = TimePoint(y, w)
                    if !haskey(tp_demand_reduction[(station, node, blk)], tp)
                        tp_demand_reduction[(station, node, blk)][tp] = 0.0
                    end
                    tp_demand_reduction[(station, node, blk)][tp] += generation
                end
            end
        end
    end

    fixed = deepcopy(demand)

    for (tp, f) in fixed
        for (node, blk) in keys(f)
            f[(node, blk)] = 0.0
        end
    end

    for ((station, node, blk), vals) in tp_demand_reduction
        for (tp, val) in vals
            f = fixed[tp]
            if haskey(f, (node, blk))
                # Note: we account for outages in the fixed stations
                if (station, blk) ∈ keys(outage_data[tp])
                    f[(node, blk)] += (val - outage_data[tp][(station, blk)]) * T[tp][blk]
                else
                    f[(node, blk)] += val * T[tp][blk]
                end
            end
        end
    end
    return fixed
end
