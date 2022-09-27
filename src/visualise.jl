#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

using Random

"""
    plot_storage(
        results::Array{Array{Dict{Symbol,Any},1},1},
        filename::String = randstring(6);
        reservoirs::Array{Symbol,1} = collect(keys(results[1][1][:reslevel])),
    )

This function plots the storage in each of the reservoirs (Mm³) alongside the total
storage (GWh) as an interactive webpage.

### Required Arguments
`results` is the output from a JADE simulation.

### Keyword Arguments
`filename` is the name of the file in which you wish to store the visualisation; the filename will be a random string if omitted.

`reservoirs` is an array of Symbols for the reservoirs that are to be visualised.
"""
function plot_storage(
    results::Array{Array{Dict{Symbol,Any},1},1},
    filename::String = randstring(6);
    reservoirs::Any = collect(keys(results[1][1][:reslevel])),
    fuels::Any = collect(keys(results[1][1][:fuel_storage])),
)
    plt = SDDP.SpaghettiPlot(results)
    for r in reservoirs
        if typeof(r) <: JuMP.Containers.DenseAxisArrayKey
            name = string(r[1])
        else
            name = string(r)
        end
        SDDP.add_spaghetti(plt; title = name) do data
            return data[:reslevel][r].out
        end
    end
    SDDP.add_spaghetti(plt; title = "Total Storage", ymax = 4.5 * 10^3) do data
        return data[:total_storage]
    end
    for f in fuels
        if typeof(f) <: JuMP.Containers.DenseAxisArrayKey
            name = string(f[1])
        else
            name = string(f)
        end
        SDDP.add_spaghetti(plt; title = name) do data
            return data[:fuel_storage][f].out
        end
    end

    SDDP.add_spaghetti(plt; title = "Historical year") do data
        return data[:inflow_year]
    end
    return SDDP.plot(plt, "storage_" * filename * ".html", open = true)
end

"""
    plot_prices(
        results::Array{Array{Dict{Symbol,Any},1},1},
        filename::String = randstring(6);
        nodes::Array{Symbol,1} = results[1][1][:prices].axes[1],
        blocks::Array{Symbol,1} = results[1][1][:prices].axes[2],
    )

This function plots the prices for each of the nodes and load blocks alongside the
total storage (GWh) as an interactive webpage.

### Required Arguments
`results` is the output from a JADE simulation.

### Keyword Arguments
`nodes` is an array of Symbols for the nodes that are to be visualised.

`filename` is the name of the file in which you wish to store the visualisation; the filename will be a random string if omitted.

`blocks` is an array of Symbols for the load blocks that are to be visualised.
"""
function plot_prices(
    results::Array{Array{Dict{Symbol,Any},1},1},
    filename::String = randstring(6);
    nodes::Array{Symbol,1} = results[1][1][:prices].axes[1],
    blocks::Array{Symbol,1} = results[1][1][:prices].axes[2],
)
    plt = SDDP.SpaghettiPlot(results)
    for n in nodes
        for b in blocks
            SDDP.add_spaghetti(plt; title = string(n) * " - " * string(b)) do data
                return data[:prices][n, b]
            end
        end
    end
    SDDP.add_spaghetti(plt; title = "Total Storage", ymax = 4.5 * 10^3) do data
        return data[:total_storage]
    end
    SDDP.add_spaghetti(plt; title = "Historical year") do data
        return data[:inflow_year]
    end
    return SDDP.plot(plt, "prices_" * filename * ".html", open = true)
end

"""
    plot_watervalues(
        results::Array{Array{Dict{Symbol,Any},1},1},
        filename::String = randstring(6);
        reservoirs::Array{Symbol,1} = collect(keys(results[1][1][:reslevel])),
    )

This function plots the water values (\$/MWh) and storage in each of the reservoirs (Mm³) alongside the total
storage (GWh) as an interactive webpage.

### Required Arguments
`results` is the output from a JADE simulation.

### Optional Arguments
`filename` is the name of the file in which you wish to store the visualisation; the filename will be a random string if omitted.

`reservoirs` is an array of Symbols for the reservoirs that are to be visualised.
"""
function plot_watervalues(
    results::Array{Array{Dict{Symbol,Any},1},1},
    filename::String = randstring(6);
    reservoirs::Any = collect(keys(results[1][1][:reslevel])),
    fuels::Any = collect(keys(results[1][1][:fuel_storage])),
)
    plt = SDDP.SpaghettiPlot(results)
    for r in reservoirs
        if typeof(r) <: JuMP.Containers.DenseAxisArrayKey
            name = string(r[1])
            rr = r[1]
        else
            name = string(r)
            rr = r
        end
        SDDP.add_spaghetti(plt; title = name * " storage") do data
            return data[:reslevel][r].out
        end
        SDDP.add_spaghetti(plt; title = name * " MWV") do data
            return data[:mwv][rr]
        end
    end
    SDDP.add_spaghetti(plt; title = "Total Storage", ymax = 4.5 * 10^3) do data
        return data[:total_storage]
    end

    for f in fuels
        if typeof(f) <: JuMP.Containers.DenseAxisArrayKey
            name = string(f[1])
            ff = f[1]
        else
            name = string(f)
            ff = f
        end
        SDDP.add_spaghetti(plt; title = name * " storage") do data
            return data[:fuel_storage][f].out
        end
        SDDP.add_spaghetti(plt; title = name * " MFV") do data
            return data[:mfv][ff]
        end
    end

    SDDP.add_spaghetti(plt; title = "Historical year") do data
        return data[:inflow_year]
    end
    return SDDP.plot(plt, "wv_" * filename * ".html", open = true)
end
