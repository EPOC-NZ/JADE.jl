#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

module TestDataValidation

using Test
using JADE

function runtests()
    jade_dir = get(ENV, "JADE_DIR", nothing)
    ENV["JADE_DIR"] = joinpath(dirname(@__DIR__), "test")
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    if jade_dir === nothing
        delete!(ENV, "JADE_DIR")
    else
        ENV["JADE_DIR"] = jade_dir
    end
    return
end

function _input_file(s)
    return joinpath(JADE.@JADE_DIR, "Input", "test1", s)
end

function test_data_thermal_stations()
    stations =
        JADE.getthermalstations(_input_file("thermal_stations.csv"), [:NI, :HAY, :SI])
    @test length(stations) == 13
    stratford = stations[:STRATFORD_PEAKERS]
    @test stratford.fuel == :GAS
    @test stratford.heatrate == 9.5
    return
end

function test_data_thermal_fuel_costs()
    costs, fuels = JADE.getfuelcosts(_input_file("thermal_fuel_costs.csv"))
    @test length(fuels) == 3
    @test costs[JADE.TimePoint(2005, 10)][:COAL] == 4.0
    @test costs[JADE.TimePoint(2009, 52)][:DIESEL] == 24.33024581
    @test length(costs) == 260
    return
end

function test_data_reservoirs()
    reservoirs = JADE.initialisereservoirs(
        _input_file("reservoirs.csv"),
        _input_file("reservoir_limits.csv"),
    )
    @test length(reservoirs) == 7
    benmore = reservoirs[:LAKE_BENMORE]
    @test benmore.initial == 322.00032233399
    return
end

function test_data_transmission_lines()
    lineoutage, arcs = JADE.gettimeseries(_input_file("transmission_outages.csv"))
    lineoutage = JADE.duplicate_data(lineoutage, [:PEAK, :SHOULDER, :OFFPEAK])

    arcs, losses = JADE.gettransarcs(
        _input_file("transmission.csv"),
        lineoutage,
        :none,
        [:PEAK, :SHOULDER, :OFFPEAK],
    )
    @test length(arcs) == 2
    if haskey(arcs, (:NI, :HAY))
        @test arcs[(:NI, :HAY)].poscapacity == 1000.0
    else
        @test arcs[(:HAY, :NI)].poscapacity == 1000.0
    end
    return
end

function test_data_hydro_arcs()
    arcs = JADE.getnaturalarcs(_input_file("hydro_arcs.csv"))
    @test length(arcs) == 17
    @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)].minflow == 0.0
    @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)].maxflow == Inf
    return
end

function test_data_hydro_stations()
    stations, station_arcs =
        JADE.gethydros(_input_file("hydro_stations.csv"), [:NI, :HAY, :SI])
    @test length(stations) == 26
    tokaanu = stations[:TOKAANU]
    @test tokaanu.node == :NI
    @test tokaanu.capacity == 240.0
    return
end

function test_data_station_outages()
    outages, stations = JADE.gettimeseries(_input_file("station_outages.csv"))
    @test length(outages) == 260
    return
end

function test_data_terminal_water_value()
    eqn = JADE.getterminalvalue(_input_file("terminal_water_value.csv"))
    @test length(eqn) == 5
    @test eqn[1].coefficient == 137.218398630355
    return
end

function test_data_demand_durations()
    duration, blocks = JADE.gettimeseries(_input_file("hours_per_block.csv"))
    demand, nodes = JADE.getdemand(_input_file("demand.csv"), duration)
    @test length(nodes) == 3
    @test setdiff(nodes, [:NI, :SI, :HAY]) == []
    @test length(demand) == 260
    @test length(keys(demand[JADE.TimePoint(2005, 1)])) == 9
    @test isapprox(demand[JADE.TimePoint(2005, 1)][(:SI, :OFFPEAK)], 60764.5, atol = 0.1)
    @test length(blocks) == 3
    @test setdiff(blocks, [:OFFPEAK, :PEAK, :SHOULDER]) == []
    @test length(duration) == 260
    @test length(keys(duration[JADE.TimePoint(2005, 1)])) == 3
    @test duration[JADE.TimePoint(2005, 1)][:OFFPEAK] == 51
    return
end

function test_data_lost_load()
    # TODO(odow): why is this commented out?
    # lost_load_cost, sectors = JADE.getshedding(_input_file("lost_load.csv"))
    # @test length(sectors) == 3
    # @test setdiff(sectors, [:INDUSTRIAL, :COMMERCIAL, :RESIDENTIAL]) == []
    # @test length(lost_load_cost[(:NI, :INDUSTRIAL)]) == 3
    return
end

end  # module

TestDataValidation.runtests()
