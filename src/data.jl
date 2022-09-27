#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

import Statistics: mean, max

"""
	parsefile(f::Function, file::String, trim::Bool = false)

This function reads in a data file and extracts each of the non-comment / non-empty rows,
splits into into a vector using ',' as the delimiter, and then applies some custom processing `f`.

### Required Arugments
`f` custom function that processes a line of the data file.

`file` path of the file.

### Keyword Arguments
`trim` if `true` trailing commas will be discarded.

"""
function parsefile(f::Function, file::String, trim::Bool = false)
    open(file, "r") do io
        while !eof(io)
            items, skip_line = getitems(io, trim)
            skip_line && continue
            f(items)
        end
    end
end

"""
	getitems(io::IOStream, trim::Bool)

This function takes an `IOStream` for a CSV file and reading the first line, decarding empty / comment rows, and
otherwise splitting into into a vector using ',' as the delimiter. This vector is returned.

### Required Arugments
`io` `IOStream` object for a CSV file.

`trim` if `true` trailing commas will be discarded.

"""
function getitems(io::IOStream, trim::Bool)
    line = strip(readline(io))
    if (line == "" || first(line) == '%') # comment char
        return String[], true
    end
    # drop trailling comments
    line = split(line, '%')[1]

    items = strip.(split(line, ","))

    if trim
        while length(items) > 0 && items[end] == ""
            deleteat!(items, length(items))
        end
    end

    if length(items) == 0 || maximum(length.(items)) == 0
        return String[], true
    end
    return items, false
end

"""
	parse_run_options(file::String)

This function loads a run file (CSV) creates a `Dict` that has keys corresponding to various
`JADE` `RunData` fields, and their associated values.

### Required Arugments
`file` is the path to the run file that is being loaded.

"""
function parse_run_options(file::String)
    run_options = Dict{String,String}()
    parsefile(file) do items
        @assert length(items) == 2 # must be [key,value]
        if haskey(run_options, items[1])
            error("Run option $(items[1]) given twice!")
        end
        return run_options[items[1]] = items[2]
    end
    return run_options
end

str2sym(s::AbstractString) = Symbol(uppercase(String(s)))

include("time.jl")

"""
    gettimeseries(file::String)

Read timeseries with arbitrary number of columns from `file`.

### Required Arguments
`file` is the name of the file from which JADE is loading Timeseries data.

# Example File

    % Data (A,B,C)
    YEAR,WEEK,A,B,C
    2005,1,4,18.68,4.49
    2005,2,4,18.68,4.49
"""
function gettimeseries(file::String)
    data = Dict{Symbol,Float64}[]
    columns = Symbol[]
    start_time = nothing
    parsefile(file, true) do items
        @assert length(items) >= 3
        if lowercase(items[1]) == "year"
            # header row
            for it in items[3:end]
                if it != ""
                    push!(columns, str2sym(it))
                end
            end
        else
            if length(data) == 0
                start_time = TimePoint(parse(Int, items[1]), parse(Int, items[2]))
            else
                current_time = TimePoint(parse(Int, items[1]), parse(Int, items[2]))
                if current_time != start_time + length(data)
                    error("Weeks in " + file + " must be contiguous")
                end
            end
            d = Dict{Symbol,Float64}()
            for (i, col) in enumerate(columns)
                d[col] = parse(Float64, items[i+2])
            end
            push!(data, d)
        end
    end
    return TimeSeries{Dict{Symbol,Float64}}(start_time, data), columns
end

#---------------------------------------------------------------------
# Power system nodes
#---------------------------------------------------------------------
struct NodeHas
    thermal::Array{Symbol}
    hydro::Array{Symbol}
end

function process_set_items(
    input::Union{String,SubString},
    setitems::Union{Vector{Symbol},Vector{Int}},
)
    function cnvt(value::Union{String,SubString})
        if typeof(setitems) == Vector{Int}
            return parse(Int, value)
        else
            return Symbol(value)
        end
    end

    if input == ""
        error("No items specifed from set " * string(setitems))
    end
    items = split(input, ";")
    for i in 1:length(items)
        items[i] = uppercase(items[i])
    end
    selected = nothing

    if typeof(setitems) == Vector{Int}
        selected = Int[]
    else
        selected = Symbol[]
    end

    for i in items
        if uppercase(i) == "ALL"
            i = string(setitems[1]) * "-" * string(setitems[end])
        end
        if i[1] != '!'
            range = split(i, "-")
            if length(range) == 1
                if cnvt(i) in setitems
                    if cnvt(i) ∉ selected
                        push!(selected, cnvt(i))
                    else
                        @warn("Duplicate value: " * i)
                    end
                else
                    error("Unknown set item " * i * " from " * string(setitems))
                end
            elseif length(range) == 2
                index1 = findfirst(isequal(cnvt(range[1])), setitems)
                index2 = findfirst(isequal(cnvt(range[2])), setitems)
                if index1 == nothing || index2 == nothing || index2 < index1
                    error("Invalid range specified: " * i)
                end
                for j in index1:index2
                    if setitems[j] ∉ selected
                        push!(selected, setitems[j])
                    else
                        @warn("Duplicate value: " * string(setitems[j]))
                    end
                end
            else
                error("Only two ends of a range can be specified: " * i)
            end
        else
            range = split(i[2:end], "-")
            if length(range) == 1
                if cnvt(i[2:end]) ∈ selected
                    deleteat!(selected, findfirst(isequal(cnvt(i[2:end])), selected))
                else
                    @warn("Attempted to exclude item that was not included: " * i[2:end])
                end
            elseif length(range) == 2
                index1 = findfirst(isequal(cnvt(range[1])), setitems)
                index2 = findfirst(isequal(cnvt(range[2])), setitems)
                if index1 == nothing || index2 == nothing || index2 <= index1
                    error("Invalid range specified: " * i)
                end
                for j in index1:index2
                    if setitems[j] ∈ selected
                        deleteat!(selected, findfirst(isequal(setitems[j]), selected))
                    else
                        @warn(
                            "Attempted to exclude item that was not included: " *
                            string(setitems[j])
                        )
                    end
                end
            else
                error("Only two ends of a range can be specified: " * i)
            end
        end
    end

    if length(selected) == 0
        error("No items selected from set " * string(setitems))
    end
    return selected
end

include(joinpath("data", "demand.jl"))
include(joinpath("data", "hydro.jl"))
include(joinpath("data", "thermal.jl"))
include(joinpath("data", "outages.jl"))
include(joinpath("data", "fixed.jl"))
include(joinpath("data", "transmission.jl"))
include(joinpath("data", "checks.jl"))
include("network.jl")

"""
    getnodes(NODES::Vector{Symbol},
             thermal_stations::Dict{Symbol,ThermalStation},
             hydro_stations::Dict{Symbol,HydroStation})

Assigns thermal and hydro plants to their nodes.

### Required Arguments
`NODES` is a vector of Symbols representing the nodes in the network.
`thermal_stations` is a vector of `ThermalStation` objects.
`hydro_stations` is a vector of `HydroStation` objects.
"""
function getnodes(
    NODES::Vector{Symbol},
    thermal_stations::Dict{Symbol,ThermalStation},
    hydro_stations::Dict{Symbol,HydroStation},
)
    nodeproperties = Dict{Symbol,NodeHas}()
    for n in NODES
        nodeproperties[n] = NodeHas(Symbol[], Symbol[])
        for (name, station) in thermal_stations
            if station.node == n
                push!(nodeproperties[n].thermal, name)
            end
        end
        for (name, station) in hydro_stations
            if station.node == n
                push!(nodeproperties[n].hydro, name)
            end
        end
    end
    return nodeproperties
end

#---------------------------------------------
# Data type definitions and functions
#--------------------------------------------
"""
An object containing all the data required to run the JADE model.

### Fields

`rundata` Run-related constants.

`thermal_stations` Dictionary of thermal power station properties.

`hydro_stations` Dictionary of hydro station properties.

`reservoirs` Dictionary of reservoir properties.

`fuel_costs` Time series of weekly fuel costs.

`carbon_content` Carbon content for each fuel type.

`inflow_mat` Structure for storing inflow scenarios.

`station_arcs` Properties of arcs that correspond to hydro releases or spills.

`natural_arcs` Arcs not going through hydro stations.

`nodehas` Dictionary of node properties.

`spMax` The maximum specific power value.

`transmission` Dictionary of transmission line properties.

`loops` Array of independent cycles in the transmission network, each entry being an array of arcs.

`durations` Timeseries with block durations.

`demand` Demand data.

`outage` Outage each week for each station.

`dr_tranches` Tranches for power-based demand response and load shedding.

`en_tranches` Tranches for energy-based demand response.

`terminal_eqns` Equations for terminal water value.

`sets` Structure containing all sets for the JADE model.
"""
mutable struct JADEData
    rundata::RunData
    thermal_stations::Dict{Symbol,ThermalStation}
    hydro_stations::Dict{Symbol,HydroStation}
    reservoirs::Dict{Symbol,Reservoir}
    fuel_costs::TimeSeries{Dict{Symbol,Float64}}
    carbon_content::Dict{Symbol,Float64}
    inflow_mat::Vector{Dict{Symbol,Vector{Float64}}}
    station_arcs::Dict{NTuple{2,Symbol},StationArc}
    natural_arcs::Dict{NTuple{2,Symbol},NaturalArc}
    nodehas::Dict{Symbol,NodeHas}
    spMax::Float64
    transmission::Dict{NTuple{2,Symbol},TransArc}
    loops::Vector{Vector{NTuple{2,Symbol}}}
    durations::TimeSeries{Dict{Symbol,Float64}}
    demand::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}
    fixed::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}
    outage::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}}
    dr_tranches::TimeSeries{Dict{Symbol,Dict{Symbol,Dict{Tuple{Symbol,Symbol},Tranche}}}}
    en_tranches::TimeSeries{Dict{Symbol,Dict{Tuple{Symbol,Vector{Symbol}},Vector{Tranche}}}}
    terminal_eqns::Array{LinearEquation}
    sets::Sets
end

"""
    JADEdata(rundata::RunData)

# Description

This function builds a `JADEdata` object, which holds all parameters for the
DOASA stage problem.

# Required Arguments

    rundata::RunData     A structure containing JADE model parameters
"""
function JADEdata(rundata::RunData)
    rundata = deepcopy(rundata)

    filedir(x) = get_file_directory(x, rundata)

    sets = Sets() # initialise JADE sets

    @info("Input load blocks") # Get block durations
    durations, sets.BLOCKS = gettimeseries(filedir("hours_per_block.csv"))

    @info("Input demand")
    demand, sets.NODES = getdemand(filedir("demand.csv"), durations)
    checkdemands(durations, demand, rundata)

    @info("Input thermal stations")
    thermal_stations = getthermalstations(filedir("thermal_stations.csv"), sets.NODES)
    sets.THERMALS = collect(keys(thermal_stations))

    @info("Input hydro stations")
    hydro_stations, station_arcs = gethydros(filedir("hydro_stations.csv"), sets.NODES)
    sets.HYDROS = collect(keys(hydro_stations))
    sets.STATION_ARCS = collect(keys(station_arcs))

    # Prepare reservoir parameters
    reservoirs =
        initialisereservoirs(filedir("reservoirs.csv"), filedir("reservoir_limits.csv"))
    sets.RESERVOIRS = collect(keys(reservoirs))

    natural_arcs = getnaturalarcs(filedir("hydro_arcs.csv"))
    sets.NATURAL_ARCS = collect(keys(natural_arcs))

    for j in sets.NATURAL_ARCS
        if j[1] ∉ [sets.JUNCTIONS; sets.RESERVOIRS]
            push!(sets.JUNCTIONS, j[1])
        end
        if j[2] ∉ [sets.JUNCTIONS; sets.RESERVOIRS]
            push!(sets.JUNCTIONS, j[2])
        end
    end
    for j in sets.STATION_ARCS
        if j[1] ∉ [sets.JUNCTIONS; sets.RESERVOIRS]
            push!(sets.JUNCTIONS, j[1])
        end
        if j[2] ∉ [sets.JUNCTIONS; sets.RESERVOIRS]
            push!(sets.JUNCTIONS, j[2])
        end
    end

    # Inform user if an arc is defined in two places. The current behaviour is that
    # it will be considered as two different arcs with the same origin and destination
    # but different properties.
    for a in sets.STATION_ARCS
        if a in sets.NATURAL_ARCS
            @warn("Arc $a appears in multiple input files.")
        end
    end

    # Get our inflows, adjusted using DIA
    adjusted_inflows, firstweekinflows = adjustinflows(filedir("inflows.csv"), rundata)
    # Archive our adjusted inflows to a file
    diatofile(adjusted_inflows, rundata.data_dir, rundata.policy_dir)
    # Get inflows into the array form we will use in the JADE model
    inflow_mat = getinflows(adjusted_inflows, firstweekinflows, rundata)
    # inflow_mat = getinflows(filedir("inflows.csv"), rundata)
    sets.CATCHMENTS = union(sets.RESERVOIRS, sets.JUNCTIONS)
    # sets.CATCHMENTS_WITH_INFLOW = collect(keys(inflow_mat[1]))
    sets.CATCHMENTS_WITH_INFLOW = collect(keys(adjusted_inflows))
    @assert length(setdiff(sets.CATCHMENTS_WITH_INFLOW, sets.CATCHMENTS)) == 0
    sets.JUNCTIONS_WITHOUT_INFLOW =
        setdiff(sets.JUNCTIONS, [sets.CATCHMENTS_WITH_INFLOW; :SEA])

    checkinflows(adjusted_inflows, firstweekinflows, rundata)

    # Set the specific power value for each reservoir
    set_reservoir_sp!(reservoirs, hydro_stations, sets, station_arcs)
    spMax = maximum([reservoirs[r].sp for r in sets.RESERVOIRS])

    # Dictionary of capacities and losses for transmission arcs
    @info("Input transmission capacities, outages and losses")

    if isfile(filedir("line_outages.csv"))
        if isfile(filedir("transmission_outages.csv"))
            @warn(
                "Reading transmission outages from 'line_outages.csv', 'transmission_outages.csv' will be ignored"
            )
        end
        lineoutage = getoutages(
            filedir("line_outages.csv"),
            collect(demand.startpoint.year:(demand.startpoint+length(demand.data)-1).year),
            collect(1:WEEKSPERYEAR),
            sets.BLOCKS,
            demand,
            durations,
        )
    elseif isfile(filedir("transmission_outages.csv"))
        @info(
            "Reading transmission outages from 'transmission_outages.csv' in DOASA compatibility mode"
        )
        lineoutage, arcs = gettimeseries(filedir("transmission_outages.csv"))
        lineoutage = duplicate_data(lineoutage, sets.BLOCKS)
    else
        error("'line_outages.csv' not found in input directory.")
    end

    transmission, rundata.losses =
        gettransarcs(filedir("transmission.csv"), lineoutage, rundata.losses, sets.BLOCKS)

    sets.TRANS_ARCS = collect(keys(transmission))

    @info("Input outages")

    if isfile(filedir("generator_outages.csv"))
        if isfile(filedir("station_outages.csv"))
            @warn(
                "Reading outages from 'generator_outages.csv', 'station_outages.csv' will be ignored"
            )
        end
        outage = getoutages(
            filedir("generator_outages.csv"),
            collect(demand.startpoint.year:(demand.startpoint+length(demand.data)-1).year),
            collect(1:WEEKSPERYEAR),
            sets.BLOCKS,
            demand,
            durations,
        )
    elseif isfile(filedir("station_outages.csv"))
        @info("Reading outages from 'station_outages.csv' in DOASA compatibility mode")
        outage, stations = gettimeseries(filedir("station_outages.csv"))
        if length(
            setdiff(union(sets.THERMALS, sets.HYDROS), collect(keys(outage.data[1]))),
        ) != 0
            error("Not all generators are listed in station_outages.csv")
        end
        outage = duplicate_data(outage, sets.BLOCKS)
    else
        error("'generator_outages.csv' not found in input directory.")
    end

    checkoutages(thermal_stations, hydro_stations, outage, sets.BLOCKS, rundata)

    @info("Compute fixed generation")

    if isfile(filedir("fixed_generation.csv"))
        if isfile(filedir("fixed_stations.csv"))
            @warn(
                "Reading fixed generation from 'fixed_generation.csv', 'fixed_stations.csv' will be ignored"
            )
        end
        fixed = getfixedgeneration(
            filedir("fixed_generation.csv"),
            collect(demand.startpoint.year:(demand.startpoint+length(demand.data)-1).year),
            collect(1:WEEKSPERYEAR),
            sets.BLOCKS,
            demand,
            outage,
            durations,
        )
    elseif isfile(filedir("fixed_stations.csv"))
        @info(
            "Reading fixed generation from 'fixed_stations.csv' in DOASA compatibility mode"
        )
        fixed = computefixed(filedir("fixed_stations.csv"), demand, outage, durations)
    else
        error("'fixed_generation.csv' not found in input directory.")
    end

    # Get loops in the network
    loops = Array{NTuple{2,Symbol}}[]

    if isfile(filedir("demand_response.csv"))
        if isfile(filedir("lost_load.csv"))
            @warn(
                "Reading demand-response tranches from 'demand_response.csv', 'lost_load.csv' will be ignored"
            )
        else
            @info("Reading demand-response tranches")
        end
        dr_tranches, en_tranches, sets.SECTORS = getdemandresponse(
            get_file_directory("demand_response.csv", rundata, verbose = false),
            demand,
            durations,
            sets.NODES,
            collect(demand.startpoint.year:(demand.startpoint+length(demand.data)-1).year),
            collect(1:WEEKSPERYEAR),
            sets.BLOCKS,
        )
    elseif isfile(filedir("lost_load.csv"))
        @info(
            "Reading demand-response tranches from 'lost_load.csv' in DOASA compatibility mode"
        )
        dr_tranches, en_tranches, sets.SECTORS = getdemandresponse(
            filedir("lost_load.csv"),
            demand,
            durations,
            sets.NODES,
            collect(demand.startpoint.year:(demand.startpoint+length(demand.data)-1).year),
            collect(1:WEEKSPERYEAR),
            sets.BLOCKS,
        )
    else
        error("'demand_response.csv' / 'lost_load.csv' not found in input directory.")
    end

    @info("Input thermal fuel properties")

    fuel_costs, carbon_content = getfuelcosts(filedir("thermal_fuel_costs.csv"))
    checkfuelcosts(fuel_costs, rundata)

    @info(
        "Recording input data files in " *
        joinpath("Output", rundata.data_dir, rundata.policy_dir, "data_files")
    )
    backup_input_files(rundata)

    return JADEData(
        rundata,
        thermal_stations,
        hydro_stations,
        reservoirs,
        fuel_costs,
        carbon_content,
        inflow_mat,
        station_arcs,
        natural_arcs,
        getnodes(sets.NODES, thermal_stations, hydro_stations),
        spMax,
        transmission,
        loops,
        durations,
        demand,
        fixed,
        outage,
        dr_tranches,
        en_tranches,
        getterminalvalue(filedir("terminal_water_value.csv")),
        sets,
    )
end

function backup_input_files(rundata::RunData)
    input_filenames = [
        "hours_per_block.csv",
        "demand.csv",
        "thermal_stations.csv",
        "hydro_stations.csv",
        "reservoirs.csv",
        "reservoir_limits.csv",
        "hydro_arcs.csv",
        "inflows.csv",
        "line_outages.csv",
        "transmission_outages.csv",
        "transmission.csv",
        "generator_outages.csv",
        "station_outages.csv",
        "fixed_generation.csv",
        "terminal_water_value.csv",
        "fixed_stations.csv",
        "demand_response.csv",
        "lost_load.csv",
        "demand_response.csv",
        "lost_load.csv",
        "thermal_fuel_costs.csv",
        "thermal_fuel_storage.csv",
        "thermal_fuel_supply.csv",
    ]

    out_path =
        joinpath(@JADE_DIR, "Output", rundata.data_dir, rundata.policy_dir, "data_files")

    if !isdir(out_path)
        mkdir(out_path)
    end

    for file in input_filenames
        path_to_file = get_file_directory(file, rundata, verbose = false)

        if isfile(path_to_file)
            cp(path_to_file, joinpath(out_path, file), force = true)
        end
    end
end
