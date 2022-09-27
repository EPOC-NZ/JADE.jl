#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#--------------------------------------------------
# Prepare reservoir-related data
#--------------------------------------------------
# Arcs in water network
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

"""
    initialisereservoirs(file::String, limits::String)

# Description

Read list of reservoirs from `file`.
Capacity and initial contents in Mm³.
The columns MUST be ordered as shown below.

# Example File

    RESERVOIR, INFLOW_REGION, CAPACITY, INI_STATE
    Lake_Benmore, SI, 423451075.96799916, 322000320.233399
"""
function initialisereservoirs(file::String, limits::String)
    reservoirs = Dict{Symbol,Reservoir}()
    parsefile(file) do items
        @assert length(items) == 3 # must be 3 columns
        if lowercase(items[1]) == "reservoir"
            return
        end
        reservoir = str2sym(items[1])
        if haskey(reservoirs, reservoir)
            error("Reservoir $(reservoir) given twice.")
        end
        return reservoirs[reservoir] = Reservoir(
            JADE.TimeSeries{Float64}(JADE.TimePoint(0, 0), []),   # capacity placeholder
            parse(Float64, items[3]),   # initial state
            0.0,                         # specific power
            JADE.TimeSeries{Vector{ContingentTranche}}(JADE.TimePoint(0, 0), []),   # contingent placeholder
            length(reservoirs) + 1, # index of reservoir, used to create compatible DOASA cut files
        )
    end

    limits, sym = gettimeseries(limits)

    for (r, res) in reservoirs
        temp = Float64[]
        for i in 1:length(limits)
            push!(temp, limits[i][Symbol(string(r) * " MAX_LEVEL")])
        end
        res.capacity = TimeSeries{Float64}(limits.startpoint, temp)

        temp = Vector{ContingentTranche}[]
        j = 1
        total = zeros(Float64, length(limits))
        while Symbol(string(r) * " MIN_" * string(j) * "_LEVEL") ∈ sym || j == 1
            if j > 1 && Symbol(string(r) * " MIN_" * string(j) * "_PENALTY") ∉ sym
                error(string(r) * " contingent storage penalty missing")
            end
            for i in 1:length(limits)
                limit = 0.0
                if Symbol(string(r) * " MIN_" * string(j) * "_LEVEL") ∈ sym
                    if j == 1
                        push!(temp, ContingentTranche[])
                        limit =
                            -limits[i][Symbol(string(r) * " MIN_" * string(j) * "_LEVEL")]
                    else
                        total[i] =
                            -limits[i][Symbol(string(r) * " MIN_" * string(j) * "_LEVEL")]
                        limit =
                            total[i] + limits[i][Symbol(
                                string(r) * " MIN_" * string(j - 1) * "_LEVEL",
                            )]
                    end

                    tranche = ContingentTranche(
                        limit,
                        limits[i][Symbol(string(r) * " MIN_" * string(j) * "_PENALTY")],
                    )
                    push!(temp[i], tranche)
                else
                    push!(temp, [ContingentTranche(0.0, 0.0)])
                end
            end
            j += 1
        end
        total .-= minimum(total)
        if maximum(total) >= 1E-8
            error("The maximum contingent storage must be constant for each reservoir.")
        end

        if length(temp) != 0
            res.contingent = TimeSeries{Vector{ContingentTranche}}(limits.startpoint, temp)
        end
    end
    return reservoirs
end

#------------------------------------------------------
# Flow arcs and hydro generator data
#------------------------------------------------------
"""
    getnaturalarcs(file::String)

# Description

Read list of hydro station data from `file`.
Canals, rivers, and other means of getting water from one place to another.
Excludes power station turbines and spillways; those are covered in the hydro_stations file.
Min/max flow in cumecs.
NA in MIN_FLOW converted to 0.0
NA in MAX_FLOW converted to 99999.0
The columns MUST be ordered as shown below.

# Example File

    ORIG,DEST,MIN_FLOW,MAX_FLOW
    Lake_Wanaka,Lake_Dunstan, NA, NA
    Lake_Hawea,Lake_Dunstan, 0.00, 99999
"""
function getnaturalarcs(file::String)
    natural_arcs = Dict{NTuple{2,Symbol},NaturalArc}()
    parsefile(file, true) do items
        if length(items) ∉ [4, 6]
            error(
                "hydro_arcs.csv should have 4 or 6 columns, " *
                string(length(items)) *
                " found",
            )
        end
        if lowercase(items[1]) == "orig"
            return
        end
        od_pair = (str2sym(items[1]), str2sym(items[2]))
        if haskey(natural_arcs, od_pair)
            error("Arc $(od_pair) given twice.")
        else
            natural_arcs[od_pair] = NaturalArc(
                (lowercase(items[3]) == "na") ? 0.0 : parse(Float64, items[3]),
                (lowercase(items[4]) == "na") ? Inf : parse(Float64, items[4]),
                (length(items) == 4 || lowercase(items[5]) == "default") ? -1.0 :
                parse(Float64, items[5]),
                (length(items) == 4 || lowercase(items[6]) == "default") ? -1.0 :
                parse(Float64, items[6]),
            )
        end
    end
    return natural_arcs
end

# Hydro generators
struct HydroStation
    node::Symbol
    capacity::Float64
    sp::Float64
    omcost::Float64
    arc::NTuple{2,Symbol}
end

"""
    gethydros(file::String)

# Description

Read list of hydro station data from `file`.
Capacity in MW;
Specific power in MW/cumec;
Spillway max flow in cumecs (can be zero if no spillway exists, or "na" for unlimited spill).
The columns MUST be ordered as shown below.

# Example File

    GENERATOR,HEAD_WATER_FROM,TAIL_WATER_TO,POWER_SYSTEM_NODE,CAPACITY,SPECIFIC_POWER,SPILLWAY_MAX_FLOW
    Arapuni,Lake_Arapuni,Lake_Karapiro,NI,196.7,0.439847649,99999
"""
function gethydros(file::String, nodes::Vector{Symbol})
    hydros = Dict{Symbol,HydroStation}()
    station_arcs = Dict{NTuple{2,Symbol},StationArc}()

    mode = -1

    parsefile(file, true) do items
        if length(items) ∉ [7, 8, 9]
            error(
                "hydro_stations.csv should have 7, 8, 9 columns, " *
                string(length(items)) *
                " found",
            )
        end

        if lowercase(items[1]) == "generator"
            if length(items) == 7
                mode = 0
            elseif length(items) == 8
                if lowercase(items[8]) == "om_cost"
                    mode = 1
                elseif lowercase(items[8]) == "overspill_penalty"
                    mode = 2
                end
            elseif length(items) == 9 &&
                   lowercase(items[8]) == "om_cost" &&
                   lowercase(items[9]) == "overspill_penalty"
                mode = 3
            end
            if mode == -1
                error(
                    "Invalid column headings, if there are 8 or 9 columns the final two must be labelled 'om_cost' and 'overspill_penalty'",
                )
            end
            return
        end
        generator = str2sym(items[1])
        if haskey(hydros, generator)
            error("Generator $(generator) given twice.")
        end
        station_arc = (str2sym(items[2]), str2sym(items[3]))
        if str2sym(items[4]) in nodes
            hydros[generator] = HydroStation(
                str2sym(items[4]),
                parse(Float64, items[5]),
                parse(Float64, items[6]),
                (mode == 0 || mode == 2) ? 0.0 : parse(Float64, items[8]),
                station_arc,
            )
        else
            error("Node " * items[4] * " for generator " * items[1] * " not found")
        end

        if haskey(station_arcs, station_arc)
            error("Station arc $(station_arc) already given.")
        else
            spillwayflow =
                station_arcs[station_arc] = StationArc(
                    (lowercase(items[7]) == "na") ? Inf : parse(Float64, items[7]),
                    (mode == 0 || mode == 1) ? -1.0 :
                    ((mode == 2) ? parse(Float64, items[8]) : parse(Float64, items[9])),
                    generator,
                )
        end
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
            mean(inflow for (tp, inflow) in zip(tps, infl) if tp.week == t) for t in weeks
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
            mean(W[i] for (i, tp) in enumerate(tps) if tp.week == t) for t in weeks
        ]
        d = [  # mean + the bigger deviation
            max(0.0, α[tp.week] + (W[i] - m[tp.week]) / √ω) for (i, tp) in enumerate(tps)
        ]
        d_mean = [  # mean of adjusted values
            mean(d[i] for (i, tp) in enumerate(tps) if tp.week == t) for t in weeks
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
    outpath = joinpath(@JADE_DIR, "Output", outpath)
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
