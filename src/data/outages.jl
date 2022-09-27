#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

function getoutages(
    outage_file::String,
    years::Vector{Int},
    weeks::Vector{Int},
    loadblocks::Vector{Symbol},
    demand::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    T::TimeSeries{Dict{Symbol,Float64}},
)
    outages = Dict{Tuple{Symbol,Symbol},Float64}[]

    for k in keys(demand)
        push!(outages, Dict{Tuple{Symbol,Symbol},Float64}())
    end

    outages = TimeSeries(demand.startpoint, outages)

    parsefile(outage_file, true) do items
        @assert length(items) == 5
        if lowercase(items[1]) == "year"
            return
        end

        years_ = process_set_items(items[1], years)
        weeks_ = process_set_items(items[2], weeks)
        loadblocks_ = process_set_items(items[3], loadblocks)
        name = str2sym(items[4])
        outage = parse(Float64, items[5])

        for blk in loadblocks_
            for y in years_
                for w in weeks_
                    tp = TimePoint(y, w)
                    if haskey(outages[tp], (name, blk))
                        error(
                            "Outage for " *
                            string((name, blk)) *
                            " in week " *
                            string(tp) *
                            " given multiple times",
                        )
                    end
                    outages[tp][(name, blk)] = outage
                end
            end
        end
    end

    return outages
end
