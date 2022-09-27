#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

import LinearAlgebra: nullspace

"""
	out_neighbors(vertex::Symbol, edges::Vector{NTuple{2,Symbol}})

Helper function to work out nodes that are adjacent and downstream.

### Required Arguments
`vertex` is the vertex in the network we are finding downstream neighbors for.

`edges` is a vector of all the edges in the hydro network.
"""
function out_neighbors(vertex::Symbol, edges::Vector{NTuple{2,Symbol}})
    neighbors = Symbol[]
    for e in edges
        if e[1] == vertex
            push!(neighbors, e[2])
        end
    end
    return neighbors
end

"""
    hasdownstream(
        sets::Sets,
        station_arcs::Dict{NTuple{2,Symbol},StationArc},
    )

This function is used to determine which hydro stations are present downstream
from all reservoirs.

### Inputs
  sets                      A JADE set structure. Should include all arcs and nodes.
  station_arcs              A dictionary we will use to get the station name from a station arc.

### Returns
  reservoir_has_downstream  A dictionary indexed by reservoirs that stores a
                            list of hydro stations that each reservoir has
                            between itself and the sea.
"""
function hasdownstream(sets::Sets, station_arcs::Dict{NTuple{2,Symbol},StationArc})
    reverseflow = Symbol[]
    @assert !isempty(sets.RESERVOIRS)
    @assert !isempty(sets.STATION_ARCS)

    # All arcs put together
    arcs = union(sets.NATURAL_ARCS, sets.STATION_ARCS)

    station_arcs2 = Dict{NTuple{2,Symbol},StationArc}()

    for (pair, arc) in station_arcs
        if arc.station ∉ reverseflow
            station_arcs2[pair] = arc
        end
    end

    for (pair, arc) in station_arcs
        if arc.station in reverseflow
            deleteat!(arcs, findall(x -> x == pair, arcs))
        end
    end

    # Structure to store neighbours of each node
    neighbors = Dict{Symbol,Vector{Symbol}}()
    for c in sets.CATCHMENTS
        neighbors[c] = out_neighbors(c, arcs)
    end

    catchments = Symbol[]
    reservoirs = Symbol[]
    for c in sets.CATCHMENTS
        if length(neighbors[c]) != 0 || c == :SEA
            push!(catchments, c)
        end
    end

    for r in sets.RESERVOIRS
        if (r in sets.CATCHMENTS && length(neighbors[r]) != 0) || r ∉ sets.CATCHMENTS
            push!(reservoirs, r)
        end
    end

    # Dictionary structure to store which stations each reservoir has downstream
    reservoir_has_downstream = Dict{Symbol,Array{Symbol}}()

    # For every reservoir, we do a depth-first search to look for hydro-stations downstream.
    for r in reservoirs

        # Colour of vertices to track progress when we traverse: 0 = white, 1 = grey, 2 = black
        node_col = Dict{Symbol,Int}()

        for n in catchments
            # Flag everything as white, unvisited
            push!(node_col, n => 0)
        end

        stop = false
        current = r
        node_col[current] = 1
        # Temporary list of hydro-stations downstream
        theList = Symbol[]
        # List of visited vertices that have a grey colouring
        grey = [r]

        # Start a depth first search
        while !stop

            # Stop going down if we reached sea. If not everything goes to sea, need to use if neighbors[current] == [] instead
            if current == :SEA
                node_col[current] = 2
                pop!(grey)
                # Move down
                current = grey[end]

            else
                # If there is anything below, loop over those neighbours
                for n in neighbors[current]

                    # The arc we might move on
                    thearc = (current, n)
                    # If it's a station arc and we haven't added the station, add it
                    if haskey(station_arcs2, thearc)
                        station = station_arcs2[thearc].station
                        if !(station in theList) && !(station in reverseflow)
                            push!(theList, station)
                        end
                    end

                    # Go to the first node unvisited so far
                    if node_col[n] == 0

                        # Move down a level
                        current = n
                        node_col[current] = 1
                        push!(grey, current)

                        break

                        # Move up if everything below has been visited
                    elseif (node_col[n] == 2) && (n == neighbors[current][end])
                        # Move up
                        if current == r
                            stop = true
                            reservoir_has_downstream[r] = theList
                        else
                            node_col[current] = 2
                            pop!(grey)
                            current = grey[end]
                        end
                    end
                end # for
            end # if
        end
    end # next reservoir

    return reservoir_has_downstream
end
