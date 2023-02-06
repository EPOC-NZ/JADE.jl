#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

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
    getthermalstations(filename::String, nodes::Vector{Symbol})

# Description

Read thermal station list from `filename`.

# Example File

```raw
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
function getthermalstations(filename::String, nodes::Vector{Symbol})
    thermal_stations = Dict{Symbol,ThermalStation}()
    for row in CSV.Rows(
        filename;
        missingstring = ["NA", "na", "default"],
        stripwhitespace = true,
        comment = "%",
    )
        row = _validate_and_strip_trailing_comment(
            row,
            [
                :GENERATOR,
                :NODE,
                :FUEL,
                :HEAT_RATE,
                :CAPACITY,
                :OMCOST,
                :START_YEAR,
                :START_WEEK,
                :END_YEAR,
                :END_WEEK,
            ],
        )
        station = str2sym(row.GENERATOR)
        if haskey(thermal_stations, station)
            error("Thermal Station ($station) already given.")
        end
        node = str2sym(row.NODE)
        if !(node in nodes)
            error("Node $node for generator $station not found")
        end
        thermal_stations[station] = ThermalStation(
            node,
            str2sym(row.FUEL),
            parse(Float64, row.HEAT_RATE),
            parse(Float64, row.CAPACITY),
            parse(Float64, row.OMCOST),
            TimePoint(parse(Int, row.START_YEAR), parse(Int, row.START_WEEK)),
            TimePoint(parse(Int, row.END_YEAR), parse(Int, row.END_WEEK)),
        )
    end
    return thermal_stations
end

"""
    getfuelcosts(file::String)

# Description

Read costs and carbon content for fuels of thermal plant.

# Example File

    % Fuel costs (\$/GJ except CO2 in \$/tonne) and carbon content (tonnes CO2/GJ)
    ,,coal,diesel,gas,CO2
    CO2_CONTENT,,0.0912,0,0.0528,
    YEAR,WEEK,,,,
    2008,1,4,33.11,5.57,0
    2008,2,4,33.11,5.57,0
"""
function getfuelcosts(filename::String)
    start_time, data = nothing, Dict{Symbol,Float64}[]
    rows = CSV.Rows(
        filename;
        missingstring = ["NA", "na", "default"],
        stripwhitespace = true,
        comment = "%",
    )
    row, row_state = iterate(rows)
    fuels = Dict(
        str2sym("$k") => parse(Float64, row[k]) for
        k in CSV.getnames(row) if !(k in (:Column1, :Column2, :CO2))
    )
    # Skip YEAR,WEEK,... row
    _, row_state = iterate(rows, row_state)
    while (ret = iterate(rows, row_state)) !== nothing
        row, row_state = ret
        time = TimePoint(parse(Int, row.Column1), parse(Int, row.Column2))
        if isempty(data)
            start_time = time
        elseif time != start_time + length(data)
            error("Weeks in $filename must be contiguous")
        end
        d = Dict{Symbol,Float64}(
            str2sym("$k") => parse(Float64, row[k]) for
            k in CSV.getnames(row) if !(k in (:Column1, :Column2))
        )
        push!(data, d)
    end
    return TimeSeries{Dict{Symbol,Float64}}(start_time, data), fuels
end
