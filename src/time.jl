#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

struct TimePoint
    year::Int
    week::Int
end

function Base.:(<)(t1::TimePoint, t2::TimePoint)
    return t1.week + WEEKSPERYEAR * t1.year < t2.week + WEEKSPERYEAR * t2.year
end # (2007, 60) > (2008, 2)

Base.:(>)(t1::TimePoint, t2::TimePoint) = <(t2, t1)

function Base.:(<=)(t1::TimePoint, t2::TimePoint)
    return t1.week + WEEKSPERYEAR * t1.year <= t2.week + WEEKSPERYEAR * t2.year
end

Base.:>=(t1::TimePoint, t2::TimePoint) = <=(t2, t1)

function Base.:(==)(t1::TimePoint, t2::TimePoint)
    return t1.year == t2.year && t1.week == t2.week
end

Base.isless(t1::TimePoint, t2::TimePoint) = t1 < t2

function Base.:(-)(t1::TimePoint, t2::TimePoint)
    years = t1.year - t2.year
    weeks = t1.week - t2.week
    if weeks < 0
        years -= 1
        weeks += WEEKSPERYEAR
    end
    return WEEKSPERYEAR * years + weeks
end

# Addition and subtraction with Int
function Base.:(+)(t1::TimePoint, t2::Int)
    nwks = t1.year * WEEKSPERYEAR + t1.week + t2
    return TimePoint(div(nwks - 1, WEEKSPERYEAR), (nwks - 1) % WEEKSPERYEAR + 1)
end

# A helper function to calculate the number of weeks between timepoints
Base.:(-)(t1::TimePoint, t2::Int) = t1 + (-t2)

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

function Base.isvalid(TS::TimeSeries, t::TimePoint)
    return 0 < t - TS.startpoint + 1 <= length(TS.data)
end

# Overload the [] accessor so we can access timeseries by a timepoint
function Base.getindex(TS::TimeSeries, t::TimePoint)
    if !isvalid(TS, t)
        throw(KeyError(t))#("$t out of range.")
    end
    return getindex(TS, t - TS.startpoint + 1)
end

# or just by the index
Base.getindex(TS::TimeSeries, t::Int) = TS.data[t]

struct TimePointIterator
    startpoint::TimePoint
    endpoint::TimePoint
end

Base.length(tpi::TimePointIterator) = tpi.endpoint - tpi.startpoint + 1

function Base.iterate(tpi::TimePointIterator, state = tpi.startpoint)
    if state > tpi.endpoint
        return nothing
    end
    return (state, state + 1)
end

function Base.keys(TS::TimeSeries)
    return TimePointIterator(TS.startpoint, TS.startpoint + length(TS.data) - 1)
end
