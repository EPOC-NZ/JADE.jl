#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

struct NaturalArc
    minflow::Float64    # flow lower bound
    maxflow::Float64    # flow upper bound
    lb_penalty::Float64
    ub_penalty::Float64
end

struct StationArc
    maxflow::Float64    # spillway upper bound
    penalty::Float64
    station::Symbol     # name of the generator
end

struct ContingentTranche
    level::Float64
    penalty::Float64
end

mutable struct Reservoir
    capacity::TimeSeries{Float64}
    initial::Float64
    sp::Float64
    contingent::TimeSeries{Vector{ContingentTranche}}
    index::Int
end

function _strip_trailing_comment(x::AbstractString)
    index = findfirst('%', x)
    if index === nothing
        return x
    end
    field = strip(x[1:index-1])
    if field == "na" || field == "NA" || field == "default"
        return missing
    end
    return field
end

_strip_trailing_comment(::Missing) = missing

function _validate_and_strip_trailing_comment(row, required, optional = Symbol[])
    row_names = CSV.getnames(row)
    @assert length(required) <= length(row_names) <= length(required) + length(optional)
    for n in required
        @assert n in row_names
    end
    for n in row_names
        @assert n in required || n in optional
    end
    dict = Dict(
        name => _strip_trailing_comment(getproperty(row, name)) for
        name in CSV.getnames(row)
    )
    return (; dict...)
end

"""
    initialisereservoirs(
        reservoirs_filename::String,
        reservoir_limits_filename::String,
    )

## Description

Read list of reservoirs from `reservoirs_filename`.

Capacity and initial contents are in units of Mm^3.

## Example `reservoirs_filename`

```raw
RESERVOIR, INFLOW_REGION, CAPACITY, INI_STATE
Lake_Benmore, SI, 423451075.96799916, 322000320.233399
```

## Example `reserovoir_limits_filename`

The `MIN_(j)_LEVEL` and `MAX_(j)_PENALTY` columns are optional.

```raw
YEAR,WEEK,Lake_Hawea MAX_LEVEL,Lake_Tekapo MAX_LEVEL,Lake_Tekapo MIN_1_LEVEL,Lake_Tekapo MIN_1_PENALTY,Lake_Tekapo MIN_2_LEVEL,Lake_Tekapo MIN_2_PENALTY
2010,1,1141.95,514.1,-100,500,-191,5000
```
"""
function initialisereservoirs(
    reservoirs_filename::String,
    reservoir_limits_filename::String,
)
    reservoirs = Dict{Symbol,Reservoir}()
    for row in CSV.Rows(
        reservoirs_filename;
        missingstring = "NA",
        stripwhitespace = true,
        comment = "%",
    )
        # TODO(odow): is CAPACITY optional?
        row = _validate_and_strip_trailing_comment(
            row,
            [:RESERVOIR, :INFLOW_REGION, :INI_STATE],
            [:CAPACITY],
        )
        reservoir = str2sym(row.RESERVOIR)
        if haskey(reservoirs, reservoir)
            error("Reservoir $(reservoir) given twice.")
        end
        reservoirs[reservoir] = Reservoir(
            TimeSeries{Float64}(TimePoint(0, 0), []),  # placeholder
            parse(Float64, row.INI_STATE),
            0.0,  # specific power
            TimeSeries{Vector{ContingentTranche}}(TimePoint(0, 0), []),  # contingent placeholder
            length(reservoirs) + 1, # index of reservoir, used to create compatible DOASA cut files
        )
    end
    limits, reservoir_limit_column_names = gettimeseries(reservoir_limits_filename)
    for (r, res) in reservoirs
        res.capacity = TimeSeries{Float64}(
            limits.startpoint,
            Float64[limits[i][Symbol("$r MAX_LEVEL")] for i in 1:length(limits)],
        )
        tranches = Vector{ContingentTranche}[ContingentTranche[] for i in 1:length(limits)]
        j = 1
        total = zeros(Float64, length(limits))
        if !(Symbol("$r MIN_$(j)_LEVEL") in reservoir_limit_column_names)
            for i in 1:length(limits)
                push!(tranches[i], ContingentTranche(0.0, 0.0))
            end
        end
        while Symbol("$r MIN_$(j)_LEVEL") in reservoir_limit_column_names
            if !(Symbol("$r MIN_$(j)_PENALTY") in reservoir_limit_column_names)
                error("$r contingent storage missing MIN_$(j)_PENALTY")
            end
            for i in 1:length(limits)
                total[i] = -limits[i][Symbol("$r MIN_$(j)_LEVEL")]
                if j > 1
                    total[i] += limits[i][Symbol("$r MIN_$(j-1)_LEVEL")]
                end
                penalty = limits[i][Symbol("$r MIN_$(j)_PENALTY")]
                push!(tranches[i], ContingentTranche(total[i], penalty))
            end
            j += 1
        end
        if maximum(total) - minimum(total) >= 1e-8
            error("The maximum contingent storage is not constant for reservoir $r.")
        end
        res.contingent = TimeSeries{Vector{ContingentTranche}}(limits.startpoint, tranches)
    end
    return reservoirs
end

"""
    getnaturalarcs(file::String)

# Description

Read list of hydro station data from `file`.

Canals, rivers, and other means of getting water from one place to another.

Excludes power station turbines and spillways; those are covered in the hydro_stations file.

`MIN_FLOW` and `MAX_FLOW` are in cumecs.

NA in MIN_FLOW converted to 0.0
NA in MAX_FLOW converted to 99999.0

# Example File

```raw
ORIG,DEST,MIN_FLOW,MAX_FLOW
Lake_Wanaka,Lake_Dunstan,NA,NA
Lake_Hawea,Lake_Dunstan,0.00,NA
```

```raw
ORIG,DEST,MIN_FLOW,MAX_FLOW,LB_PENALTY,UB_PENALTY
Lake_Wanaka,Lake_Dunstan,NA,NA,100,100
Lake_Hawea,Lake_Dunstan,0.00,NA,50,50
```
"""
function getnaturalarcs(filename::String)
    natural_arcs = Dict{NTuple{2,Symbol},NaturalArc}()
    for row in CSV.Rows(
        filename;
        missingstring = ["NA", "na", "default"],
        stripwhitespace = true,
        comment = "%",
    )
        row = _validate_and_strip_trailing_comment(
            row,
            [:ORIG, :DEST, :MIN_FLOW, :MAX_FLOW],
            [:LB_PENALTY, :UB_PENALTY],
        )
        od_pair = (str2sym(row.ORIG), str2sym(row.DEST))
        if haskey(natural_arcs, od_pair)
            error("Arc $(od_pair) given twice.")
        end
        natural_arcs[od_pair] = NaturalArc(
            parse(Float64, coalesce(row.MIN_FLOW, "0")),
            parse(Float64, coalesce(row.MAX_FLOW, "Inf")),
            parse(Float64, coalesce(get(row, :LB_PENALTY, missing), "-1")),
            parse(Float64, coalesce(get(row, :UB_PENALTY, missing), "-1")),
        )
    end
    return natural_arcs
end

struct HydroStation
    node::Symbol
    capacity::Float64
    sp::Float64
    omcost::Float64
    arc::NTuple{2,Symbol}
end

"""
    gethydros(file::String)

## Description

Read list of hydro station data from `file`.
Capacity in MW;
Specific power in MW/cumec;
Spillway max flow in cumecs (can be zero if no spillway exists, or "na" for unlimited spill).

## Example File

```raw
GENERATOR,HEAD_WATER_FROM,TAIL_WATER_TO,POWER_SYSTEM_NODE,CAPACITY,SPECIFIC_POWER,SPILLWAY_MAX_FLOW
Arapuni,Lake_Arapuni,Lake_Karapiro,NI,196.7,0.439847649,na
```

```raw
GENERATOR,HEAD_WATER_FROM,TAIL_WATER_TO,POWER_SYSTEM_NODE,CAPACITY,SPECIFIC_POWER,SPILLWAY_MAX_FLOW,OM_COST,OVERSPILL_PENALTY
Arapuni,Lake_Arapuni,Lake_Karapiro,NI,196.7,0.439847649,na,100,100
```
"""
function gethydros(filename::String, nodes::Vector{Symbol})
    hydros = Dict{Symbol,HydroStation}()
    station_arcs = Dict{NTuple{2,Symbol},StationArc}()
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
                :HEAD_WATER_FROM,
                :TAIL_WATER_TO,
                :POWER_SYSTEM_NODE,
                :CAPACITY,
                :SPECIFIC_POWER,
                :SPILLWAY_MAX_FLOW,
            ],
            [:OM_COST, :OVERSPILL_PENALTY],
        )
        generator = str2sym(row.GENERATOR)
        if haskey(hydros, generator)
            error("Generator $(generator) given twice.")
        end
        station_arc = (str2sym(row.HEAD_WATER_FROM), str2sym(row.TAIL_WATER_TO))
        if haskey(station_arcs, station_arc)
            error("Station arc $(station_arc) already given.")
        end

        power_system_node = str2sym(row.POWER_SYSTEM_NODE)
        if !(power_system_node in nodes)
            error("Node $power_system_node for generator $generator not found")
        end
        hydros[generator] = HydroStation(
            power_system_node,
            parse(Float64, row.CAPACITY),
            parse(Float64, row.SPECIFIC_POWER),
            parse(Float64, get(row, :OM_COST, "0.0")),
            station_arc,
        )
        station_arcs[station_arc] = StationArc(
            parse(Float64, coalesce(row.SPILLWAY_MAX_FLOW, "Inf")),
            parse(Float64, get(row, :OVERSPILL_PENALTY, "-1")),
            generator,
        )
    end
    return hydros, station_arcs
end

#---------------------------------------------------
# Prepare inflow-related data
#---------------------------------------------------
"""
    adjustinflows(inflow_file::String, rundata::RunData)

# Description

This method applies DIA to the inflow file.
"""
function adjustinflows(inflow_file::String, rundata::RunData)
    if rundata.dialength >= 52
        error("Correlation length too large.")
    end
    # correlation length ω
    ω = if rundata.dialength <= 0
        1
    else
        rundata.dialength
    end

    catchments = Symbol[]
    regions = Symbol[]
    tmp_inflows = Dict{Symbol,Dict{TimePoint,Float64}}()
    # wrangle data
    parsefile(inflow_file, true) do items
        if lowercase(items[1]) == "catchment"
            for it in items[3:end]
                push!(catchments, str2sym(it))
                tmp_inflows[str2sym(it)] = Dict{TimePoint,Float64}()
            end
        elseif lowercase(items[1]) == "inflow_region"
            for it in items[3:end]
                push!(regions, str2sym(it))
            end
            @assert length(regions) == length(catchments)
        elseif lowercase(items[1]) == "year"
            return
        else
            year = parse(Int, items[1])
            week = parse(Int, items[2])
            tp = TimePoint(year, week)
            for (i, catchment) in enumerate(catchments)
                tmp_inflows[catchment][tp] = parse(Float64, items[i+2])
            end
        end
    end
    inflows = Dict{Symbol,TimeSeries{Float64}}()
    first_inflows = Dict{Symbol,Float64}()

    for (location, val) in tmp_inflows
        tps, infl = collect(keys(val)), collect(values(val))
        # sort all our data
        p = sortperm(tps)
        permute!(tps, p)
        permute!(infl, p)

        first_inflows[location] = val[TimePoint(rundata.start_yr, rundata.start_wk)]

        weeks = sort(unique([tp.week for tp in tps]))
        years = sort(unique([tp.year for tp in tps]))
        N = length(years)
        # Follow the steps in DOASA paper
        α = [  # mean historical inflow
            Statistics.mean(inflow for (tp, inflow) in zip(tps, infl) if tp.week == t)
            for t in weeks
        ]
        W = [] # rolling total inflow starting from week ty
        for ty in 1:length(infl)
            last = ((ty + ω - 2) % length(infl) + 1)
            if last < ty
                push!(W, sum(infl[ty:length(infl)]) + sum(infl[1:last]))
            else
                push!(W, sum(infl[ty:last]))
            end
        end
        m = [  # mean of the rolling inflow
            Statistics.mean(W[i] for (i, tp) in enumerate(tps) if tp.week == t) for
            t in weeks
        ]
        d = [  # mean + the bigger deviation
            max(0.0, α[tp.week] + (W[i] - m[tp.week]) / √ω) for (i, tp) in enumerate(tps)
        ]
        d_mean = [  # mean of adjusted values
            Statistics.mean(d[i] for (i, tp) in enumerate(tps) if tp.week == t) for
            t in weeks
        ]
        k = [  # adjusted inflows
            isnan(d[i] * α[tp.week] / d_mean[tp.week]) ? α[tp.week] :
            d[i] * α[tp.week] / d_mean[tp.week] for (i, tp) in enumerate(tps)
        ]
        inflows[location] = TimeSeries(minimum(tps), k)
    end
    return inflows, first_inflows
end

"""
    diatofile(
        adjusted::Dict{Symbol,TimeSeries{Float64}},
        outpath::String,
        policy_dir::String,
    )

# Description

This function saves DIA adjusted inflows to file.

# Required Arguments
  `adjusted` is the dictionary of DIA-adjusted inflows
  `outpath` is the filename for adjusted inflows
  `policy_dir` is the subdirectory where the outputs are stored
"""
function diatofile(
    adjusted::Dict{Symbol,TimeSeries{Float64}},
    outpath::String,
    policy_dir::String,
)
    outpath = joinpath(@__JADE_DIR__, "Output", outpath)
    if !ispath(joinpath(outpath, policy_dir))
        mkpath(joinpath(outpath, policy_dir))
    end
    open(joinpath(outpath, policy_dir * "/adjusted_inflows.csv"), "w") do io
        # print header
        print(io, "YEAR, WEEK")
        for (name, location) in adjusted
            print(io, ", $(name)")
        end
        print(io, "\n")
        # print data
        l = first(keys(adjusted))
        starttime = adjusted[l].startpoint
        endtime = starttime + length(adjusted[l]) - 1
        for t in starttime:endtime
            print(io, t.year, ", ", t.week)
            for (name, location) in adjusted
                print(io, ", ", location[t])
            end
            print(io, "\n")
        end
    end
end

function getinflows(
    adjusted::Dict{Symbol,TimeSeries{Float64}},
    firstweekinflows::Dict{Symbol,Float64},
    rundata::RunData,
)
    # A dictionary of inflow scenarios indexed by
    #   the week of year and catchment
    inflows = Dict{Symbol,Vector{Float64}}[]
    # Allocate storage
    for week in 1:WEEKSPERYEAR
        push!(inflows, Dict{Symbol,Vector{Float64}}())
        for l in keys(adjusted)
            inflows[week][l] = Float64[]
        end
    end

    for year in rundata.sample_years
        samplestart = TimePoint(year, 1)
        samplestop = TimePoint(year, WEEKSPERYEAR)

        # Push data (in order)
        for (name, location) in adjusted
            for t in samplestart:samplestop
                push!(inflows[t.week][name], location[t])
            end
        end
    end

    # locations =
    @assert length(inflows[1][first(keys(adjusted))]) == rundata.nscenarios

    if rundata.first_week_known
        # In start week, only consider inflows from "problem start year"
        #   (DOASA feature)
        for (name, location) in adjusted
            inflows[rundata.start_wk][name] =
                repeat([firstweekinflows[name]], outer = rundata.nscenarios)
        end
    end

    return inflows
end

"""
Returns an inflow matrix (actually a vector) for a given year
"""
function getinflows_for_historical(
    inflowsfile::String,
    rundata::RunData,
    year::Union{Int,Vector{Int}},
)
    rd = deepcopy(rundata)
    rd.dialength = 1
    if typeof(year) == Int
        rd.sample_years = [year]
        rd.nscenarios = 1
    else
        rd.sample_years = year
        rd.nscenarios = length(year)
    end
    adjusted_inflows, firstweekinflows = adjustinflows(inflowsfile, rd)
    return getinflows(adjusted_inflows, firstweekinflows, rd)
end

function getSequences(file::String)
    num_weeks = 0
    sequences = Vector{Vector{Int}}()
    parsefile(file, true) do items
        if num_weeks == 0
            num_weeks = length(items)
        elseif num_weeks != length(items)
            error("All custom inflows sequences must have the same length.")
        end
        sequence = Vector{Int}()
        for i in 1:num_weeks
            push!(sequence, parse(Int, items[i]))
        end
        return push!(sequences, sequence)
    end
    println(sequences)
    return sequences
end

#---------------------------------------------------
# Get terminal water value
#---------------------------------------------------
function set_reservoir_sp!(
    reservoirs::Dict{Symbol,Reservoir},
    hydros::Dict{Symbol,HydroStation},
    sets::Sets,
    station_arcs::Dict{NTuple{2,Symbol},StationArc},
)
    # Compute which hydro stations are downstream from each reservoir
    reservoir_has_downstream = hasdownstream(sets, station_arcs)

    # Update the specific power of each reservoir:
    #   conversion factor for m^3 -> MWh
    for (name, reservoir) in reservoirs
        if name ∈ keys(reservoir_has_downstream)
            sp = sum(hydros[station].sp for station in reservoir_has_downstream[name])
            reservoir.sp = sp / 3600
        else
            reservoir.sp = 0
        end
    end
end

struct LinearEquation
    intercept::Float64
    coefficient::Float64
end

"""
    getterminalvalue(file::String)

# Description

At the final stage terminal cost is a piecewise linear function of the national
stored energy. We define a number of linear equations.

- Value of stored energy in hydro lakes at the end of the time horizon.
- STORED_ENERGY is the cumulative stored energy in GWh.
- VALUE is the marginal value in \$/MWh. This column should be a decreasing sequence.

The columns MUST be ordered as shown below.

# Example File

    STORED_ENERGY, VALUE
    1000, 137.218398630355   % first 1000 GWh is worth ~ \$137/MWh
    1500, 85.9718321526058
"""
function getterminalvalue(file::String)
    terminal_equations = LinearEquation[]
    current_value = 0.0
    current_energy = 0.0
    parsefile(file, true) do items
        @assert length(items) == 2
        if lowercase(items[1]) == "stored_energy"
            return
        end
        energy = 1000.0 * parse(Float64, items[1]) # GWh -> MWh conversion
        value = parse(Float64, items[2])

        # equation: y = current_value + (x - energy) * value
        #   intercept   = curent_value - energy * value
        #   coefficient = value
        push!(
            terminal_equations,
            LinearEquation(current_value - current_energy * value, value),
        )
        current_value += (energy - current_energy) * value
        return current_energy = energy
    end
    return [
        LinearEquation(eqn.intercept - current_value, eqn.coefficient) for
        eqn in terminal_equations
    ]
end
