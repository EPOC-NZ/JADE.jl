#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#-----------------------------------------------------
# Parameters relating to thermal stations
#-----------------------------------------------------
struct ThermalStation
    node::Symbol
    fuel::Symbol
    heatrate::Float64
    capacity::Float64
    omcost::Float64
    commission::TimePoint
    decommission::TimePoint
end

struct FuelStorage
    initial::Float64    # initial storage (TJ)
    maximum::Float64    # maximum storage (TJ)
end

"""
    getthermalstations(file::String, nodes::Vector{Symbol})

# Description

Read thermal station list from `file`.
The columns MUST be ordered as shown below.

# Example File

```
% Thermal station data,,,,,,,,,
% Heat rate in GJ/MWh,,,,,,,,,
% Capacity in MW,,,,,,,,,
% O&M cost in \$/MWh,,,,,,,,,
GENERATOR,NODE,FUEL,HEAT_RATE,CAPACITY,OMCOST,START_YEAR,START_WEEK,END_YEAR,END_WEEK
Stratford_220KV,NI,gas,7.6,377,0,0,0,0,0
Huntly_e3p,NI,gas,7.2,403,0,2007,23,0,0
Huntly_main_g1,NI,coal,10.3,250,0,0,0,0,0
Huntly_main_g2,NI,coal,10.3,250,0,0,0,0,0
```
"""
function getthermalstations(file::String, nodes::Vector{Symbol})
    thermal_stations = Dict{Symbol,ThermalStation}()
    for row in CSV.Rows(
        file;
        missingstring = "NA",
        stripwhitespace = true,
        comment = "%",
    )
        generator, node = str2sym(row.GENERATOR), str2sym(row.NODE)
        if haskey(thermal_stations, generator)
            error("Thermal Station ($generator) already given")
        elseif !(node in nodes)
            error("Node $node for generator $generator not found")
        end
        thermal_stations[generator] = ThermalStation(
            node,
            str2sym(row.FUEL),
            parse(Float64, row.HEAT_RATE),
            parse(Float64, row.CAPACITY),
            parse(Float64, row.OMCOST),
            TimePoint(
                parse(Int, row.START_YEAR),
                parse(Int, row.START_WEEK),
            ),
            TimePoint(
                parse(Int, row.END_YEAR),
                parse(Int, row.END_WEEK),
            ),
        )
    end
    return thermal_stations
end

"""
    getfuelcosts(file::String)

# Description

Read costs and carbon content for fuels of thermal plant.

# Example File

```
% Fuel costs (\$/GJ except CO2 in \$/tonne) and carbon content (tonnes CO2/GJ)
,,coal,diesel,gas,CO2
CO2_CONTENT,,0.0912,0,0.0528,
YEAR,WEEK,,,,
2008,1,4,33.11,5.57,0
2008,2,4,33.11,5.57,0
```
"""
function getfuelcosts(file::String)
    data = Dict{Symbol,Float64}[]
    fuels = Symbol[]
    CO2 = Float64[]
    start_time = nothing
    state = 0
    parsefile(file, true) do items
        if lowercase(items[1]) == "" && state == 0
            # header row
            state = 1
            for it in items[3:end]
                push!(fuels, str2sym(it))
            end
        elseif lowercase(items[1]) == "co2_content"
            # get carbon content for fuels
            # the final 'fuel cost' column is the carbon price
            for it in items[3:length(fuels)+1]
                push!(CO2, parse(Float64, it))
            end
        elseif lowercase(items[1]) == "year" && state == 1
            # year / week header row
            state = 2
        elseif state == 2
            @assert length(items) >= 3
            # fuel price data
            if length(data) == 0
                start_time = TimePoint(parse(Int, items[1]), parse(Int, items[2]))
            else
                current_time = TimePoint(parse(Int, items[1]), parse(Int, items[2]))
                if current_time != start_time + length(data)
                    error("Weeks in " + file + " must be contiguous")
                end
            end
            d = Dict{Symbol,Float64}()
            for (i, fuel) in enumerate(fuels)
                d[fuel] = parse(Float64, items[i+2])
            end
            push!(data, d)
        end
    end
    return TimeSeries{Dict{Symbol,Float64}}(start_time, data), Dict(zip(fuels, CO2))
end
