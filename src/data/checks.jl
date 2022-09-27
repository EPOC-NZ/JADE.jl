#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

function checkfuelcosts(fuelcosts::TimeSeries{Dict{Symbol,Float64}}, rundata::RunData)
    StartTime = TimePoint(rundata.start_yr, rundata.start_wk)
    if StartTime < fuelcosts.startpoint
        error("No fuel cost data for run start week.")
    end
    if StartTime + rundata.number_of_wks > fuelcosts.startpoint + length(fuelcosts)
        error("Some weeks are outside range of fuel cost data.")
    end
    return nothing
end

# Check that demands and durations match with run file
function checkdemands(
    durations::TimeSeries{Dict{Symbol,Float64}},
    demand::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    rundata::RunData,
)
    StartTime = TimePoint(rundata.start_yr, rundata.start_wk)
    if (demand.startpoint != durations.startpoint) || length(demand) != length(durations)
        error(
            "Demand and hours_per_block files must start at same week and contain same number of weeks.",
        )
    end
    if (length(demand) < rundata.number_of_wks)
        @info("Demand number of weeks = ", length(demand))
        @info("Run length = ", rundata.number_of_wks)
        error("Number of weeks of demand is less than run length.")
    end
    if StartTime < demand.startpoint
        @info("Demand start point = ", demand.startpoint)
        @info("Run start point = ", StartTime)
        error("No demand data for run start week.")
    end
    if StartTime + rundata.number_of_wks > demand.startpoint + length(demand)
        @info("Demand end point = ", demand.startpoint + length(demand) - 1)
        @info("Run end point = ", StartTime + rundata.number_of_wks - 1)
        error("Some weeks are outside range of demand data.")
    end
    return nothing
end

# Check inflows match with run file
function checkinflows(
    inflows::Dict{Symbol,TimeSeries{Float64}},
    first_inflows::Dict{Symbol,Float64},
    rundata::RunData,
)
    SampleStartTime = TimePoint(rundata.sample_years[1], 1)
    SampleEndTime = TimePoint(rundata.sample_years[length(rundata.sample_years)], 52)
    # Check that all are available
    for (location, val) in inflows
        InflowEndPoint = inflows[location].startpoint + (length(inflows[location]) - 1)
        if (inflows[location].startpoint > SampleStartTime) ||
           (InflowEndPoint < SampleEndTime)
            error(" Inflow sample years outside range of inflows file")
        end
    end
    return nothing
end

"""
This function checks that there are no outages greater than the capacity of a station and outage data spans run length
"""
function checkoutages(
    thermals::Dict{Symbol,ThermalStation},
    hydros::Dict{Symbol,HydroStation},
    outages::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    loadblocks::Vector{Symbol},
    rundata::RunData,
)
    for (name, station) in thermals
        for b in loadblocks
            for tp in keys(outages)
                if !(
                    (name, b) ∉ keys(outages[tp]) ||
                    station.capacity >= outages[tp][(name, b)]
                )
                    error(
                        string(name) *
                        " outage exceeds capacity in block " *
                        string(b) *
                        ", " *
                        string(tp),
                    )
                end
            end
        end
    end
    for (name, station) in hydros
        for b in loadblocks
            for tp in keys(outages)
                if !(
                    (name, b) ∉ keys(outages[tp]) ||
                    station.capacity >= outages[tp][(name, b)]
                )
                    error(
                        string(name) *
                        " outage exceeds capacity in block " *
                        string(b) *
                        ", " *
                        string(tp),
                    )
                end
            end
        end
    end
    StartTime = TimePoint(rundata.start_yr, rundata.start_wk)
    if StartTime < outages.startpoint
        @info("Outages start point = ", outages.startpoint)
        @info("Run start point = ", StartTime)
        error("No outage data for run start week.")
    end
    if StartTime + rundata.number_of_wks > outages.startpoint + length(outages)
        @info("Outages end point = ", outages.startpoint + length(outages) - 1)
        @info("Run end point = ", StartTime + rundata.number_of_wks - 1)
        error("Some weeks are outside range of outages data.")
    end
    return nothing
end
