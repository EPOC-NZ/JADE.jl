#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#----------------------------------------------------
# Prepare demand-related data
#----------------------------------------------------
"""
    getdemand(file::String, durations::JADE.TimeSeries{Dict{Symbol,Float64}})

# Read the demand in MW

# Example

    % Demand (MW):
    NODE,YEAR,WEEK,peak,shoulder,offpeak
    NI,2005,1,2132.15,1616.51,1116.58
"""
function getdemand(file::String, durations::JADE.TimeSeries{Dict{Symbol,Float64}})
    # want a TimeSeries[t] = Dict(
    #   (:NI, :peak) => 1.0
    # )
    data = Dict{TimePoint,Dict{Tuple{Symbol,Symbol},Float64}}()
    demand_blocks = Symbol[]
    nodes = Symbol[]
    start_year = 0
    start_week = 0
    parsefile(file, true) do items
        @assert length(items) >= 4 # at least one demand block
        if lowercase(items[1]) == "node"
            # header row
            for it in items[4:end]
                if it != ""
                    push!(demand_blocks, str2sym(it))
                end
            end
        else
            timepoint = TimePoint(parse(Int, items[2]), parse(Int, items[3]))
            if length(data) == 0
                start_year = timepoint.year
                start_week = timepoint.week
            end
            node = str2sym(items[1])
            if !(node in nodes)
                push!(nodes, node)
            end
            if !haskey(data, timepoint)
                data[timepoint] = Dict{Tuple{Symbol,Symbol},Float64}()
            end
            for (i, block) in enumerate(demand_blocks)
                data[timepoint][(node, block)] = parse(Float64, items[i+3])
            end
        end
    end
    tps, vals = collect(keys(data)), collect(values(data))
    p = sortperm(tps)

    demands = TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}(
        TimePoint(start_year, start_week),
        vals[p],
    )

    for tp in tps
        for k in keys(demands[tp])
            demands[tp][k] *= durations[tp][k[2]]
        end
    end

    return demands, nodes
end

#---------------------------------------------
# Load shedding
#---------------------------------------------
function getdemandresponse(file::String, demand, durations, nodes, years, weeks, loadblocks)
    ptranches = Dict{Symbol,Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}}[]
    etranches = Dict{Symbol,Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}}[]

    for i in 1:length(demand.data)
        push!(ptranches, Dict{Symbol,Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}}())
        push!(etranches, Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}())
    end

    ptranches = TimeSeries(demand.startpoint, ptranches)
    etranches = TimeSeries(demand.startpoint, etranches)

    sectors = Symbol[]
    parsefile(file, true) do items
        if length(items) == 10
            if lowercase(items[1]) == "demand"
                return
            end

            sector = str2sym(items[1])

            if !(sector in sectors)
                push!(sectors, sector)
            end

            name = str2sym(items[2])
            nodes_ = process_set_items(items[3], nodes)
            years_ = process_set_items(items[4], years)
            weeks_ = process_set_items(items[5], weeks)
            loadblocks_ = process_set_items(items[6], loadblocks)
            mode = lowercase(items[7])
            type_ = lowercase(items[8])
            bound = parse(Float64, items[9])
            bidprice = parse(Float64, items[10])
        else
            @assert length(items) == 7

            if lowercase(items[1]) == "node"
                return
            end

            sector = str2sym(items[3])

            if !(sector in sectors)
                push!(sectors, sector)
            end

            name = str2sym(items[4])
            nodes_ = [str2sym(items[1])]
            years_ = years
            weeks_ = weeks
            loadblocks_ = loadblocks
            mode = "power"
            type_ = "proportional"
            bound = parse(Float64, items[5]) * parse(Float64, items[6])
            bidprice = parse(Float64, items[7])
        end
        if mode == "power"
            for y in years_
                for w in weeks_
                    t = TimePoint(y, w)
                    for n in nodes_
                        if !haskey(ptranches[t], n)
                            ptranches[t][n] =
                                Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}()
                        end
                        for l in loadblocks_
                            if !haskey(ptranches[t][n], l)
                                ptranches[t][n][l] = Dict{Tuple{Symbol,Symbol},Tranche}()
                            end
                            if haskey(ptranches[t][n][l], (sector, name))
                                error("Two demand response tranches with same name found.")
                            end
                            if type_ == "absolute"
                                ptranches[t][n][l][(sector, name)] =
                                    Tranche(bound, bidprice)
                            elseif type_ == "proportional"
                                if durations[t][l] == 0
                                    ptranches[t][n][l][(sector, name)] =
                                        Tranche(0.0, bidprice)
                                else
                                    ptranches[t][n][l][(sector, name)] = Tranche(
                                        max(
                                            0.0,
                                            bound * demand[t][(n, l)] / durations[t][l],
                                        ),
                                        bidprice,
                                    )
                                end
                            else
                                error(
                                    "Unrecognised demand tranche type " *
                                    type_ *
                                    " use 'absolute' or 'proportional'",
                                )
                            end
                        end
                    end
                end
            end
        elseif mode == "energy"
            for y in years_
                for w in weeks_
                    t = TimePoint(y, w)
                    for n in nodes_
                        if !haskey(etranches[t], n)
                            etranches[t][n] =
                                Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}()
                        end
                        if !haskey(etranches[t][n], (sector, loadblocks_))
                            etranches[t][n][(sector, loadblocks_)] = Tranche[]
                        end

                        if type_ == "absolute"
                            push!(
                                etranches[t][n][(sector, loadblocks_)],
                                Tranche(bound, bidprice),
                            )
                        elseif type_ == "proportional"
                            push!(
                                etranches[t][n][(sector, loadblocks_)],
                                Tranche(
                                    max(
                                        0.0,
                                        bound * sum(
                                            demand[t][(n, l)] * durations[t][l] for
                                            l in loadblocks_
                                        ),
                                    ) / 1000,
                                    bidprice,
                                ),
                            )
                        else
                            error(
                                "Unrecognised demand tranche type " *
                                type_ *
                                " use 'absolute' or 'proportional'",
                            )
                        end
                    end
                end
            end
        else
            error("Demand-response mode must be either 'power' or 'energy'")
        end
    end
    return ptranches, etranches, sectors
end
