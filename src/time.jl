#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

import Base: (<), (>), (<=), (>=), (+), (-), (==), getindex, length, isless, keys

struct TimePoint
    year::Int
    week::Int
end

# Functions by Oscar
function Base.:<(t1::TimePoint, t2::TimePoint)
    return (t1.week + WEEKSPERYEAR * t1.year < t2.week + WEEKSPERYEAR * t2.year)
end # (2007, 60) > (2008, 2)
Base.:>(t1::TimePoint, t2::TimePoint) = <(t2, t1)
function Base.:<=(t1::TimePoint, t2::TimePoint)
    return (t1.week + WEEKSPERYEAR * t1.year <= t2.week + WEEKSPERYEAR * t2.year)
end
Base.:>=(t1::TimePoint, t2::TimePoint) = <=(t2, t1)
function ==(t1::TimePoint, t2::TimePoint)
    return ((t1.year == t2.year) && (t1.week == t2.week))
end
# ==(t1::TimePoint, t2::TimePoint) = (<=(t2, t1)) && (>=(t2, t1))

Base.isless(t1::TimePoint, t2::TimePoint) = t1 < t2

# A helper function to calculate the number of weeks between timepoints
function -(t1::TimePoint, t2::TimePoint)
    years = t1.year - t2.year
    weeks = t1.week - t2.week
    if weeks < 0
        years -= 1
        weeks += WEEKSPERYEAR
    end
    return WEEKSPERYEAR * years + weeks
end

# Addition and subtraction with Int
function +(t1::TimePoint, t2::Int)
    nwks = t1.year * WEEKSPERYEAR + t1.week + t2
    return TimePoint(div(nwks - 1, WEEKSPERYEAR), (nwks - 1) % WEEKSPERYEAR + 1)
end
# A helper function to calculate the number of weeks between timepoints
-(t1::TimePoint, t2::Int) = t1 + (-t2)

# Time series functions
struct TimeSeries{T}
    startpoint::TimePoint
    data::Vector{T}
end
Base.length(ts::TimeSeries) = length(ts.data)

function Base.iterate(ts::TimeSeries, state = (ts.startpoint, ts.data[1]))
    element, count = state

    if element == nothing
        return nothing
    end

    if !isvalid(ts, element + 1)
        return ((element, ts[element]), (nothing, nothing))
    else
        return ((element, ts[element]), (element + 1, ts[element+1]))
    end
end

function isvalid(TS::TimeSeries, t::TimePoint)
    i = t - TS.startpoint + 1
    return (i > 0 && i <= length(TS.data))
end

# Overload the [] accessor so we can access timeseries by a timepoint
function Base.getindex(TS::TimeSeries, t::TimePoint)
    if isvalid(TS, t)
        i = t - TS.startpoint + 1
        return TS.data[i]
    else
        throw(KeyError(t))#("$t out of range.")
    end
end

# or just by the index
function Base.getindex(TS::TimeSeries, t::Int)
    return TS.data[t]
end

struct TimePointIterator
    startpoint::TimePoint
    endpoint::TimePoint
end

Base.length(tpi::TimePointIterator) = tpi.endpoint - tpi.startpoint + 1

function Base.iterate(tpi::TimePointIterator, state = tpi.startpoint)
    if state > tpi.endpoint
        return nothing
    else
        return (state, state + 1)
    end
end

function keys(TS::TimeSeries)
    return TimePointIterator(TS.startpoint, TS.startpoint + length(TS.data) - 1)
end

"""
Determines if a time point is between two other time points, inclusive
"""
function between(t1::TimePoint, t2::TimePoint, t3::TimePoint)
    return t3 <= t2 && t3 >= t1
end

"""
For a data frame containing columns named WEEK and YEAR, return only the subset
of data that falls within a relevant time range.
"""
function cleandata(data::DataFrame, t1::TimePoint, t2::TimePoint)
    if !(:WEEK in names(data) && :YEAR in names(data))
        error("Data should contain columns named WEEK and YEAR.")
    else
        relevant =
            map((x) -> between(t1, t2, x), map(TimePoint, data[:, :YEAR], data[:, :WEEK]))
        return data[relevant.==true, :]
    end
end

"""
This function returns the year and week, given a stage.
"""
function calculatetime(stage::Int, starttime::TimePoint)
    nyrs = div(starttime.week + stage - 2, WEEKSPERYEAR)    # years passed
    year = starttime.year + nyrs
    week = (starttime.week + stage - 1) - nyrs * WEEKSPERYEAR
    return year, week
end
