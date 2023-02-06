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

_input_file(s) = joinpath(JADE.@__JADE_DIR__, "Input", "test1", s)

_validation_file(s) = joinpath(@__DIR__, "data_validation_files", s)

function test_data_thermal_stations()
    stations =
        JADE.getthermalstations(_validation_file("thermal_stations_small.csv"), [:NI])
    @test length(stations) == 1
    @test stations[:STRATFORD_220KV] == JADE.ThermalStation(
        :NI,
        :GAS,
        7.6,
        377.0,
        0.0,
        JADE.TimePoint(2008, 1),
        JADE.TimePoint(2019, 21),
    )
    return
end

function test_data_thermal_stations_node_error()
    filename = _validation_file("thermal_stations_small.csv")
    @test_throws(
        ErrorException("Node NI for generator STRATFORD_220KV not found"),
        JADE.getthermalstations(filename, [:SI]),
    )
    return
end

function test_data_thermal_stations_duplicate_error()
    filename = _validation_file("thermal_stations_small_duplicate.csv")
    @test_throws(
        ErrorException("Thermal Station (STRATFORD_220KV) already given."),
        JADE.getthermalstations(filename, [:NI]),
    )
    return
end

function test_data_thermal_fuel_costs()
    costs, fuels = JADE.getfuelcosts(_validation_file("thermal_fuel_costs_small.csv"))
    @test length(fuels) == 3
    @test costs[JADE.TimePoint(2019, 1)][:COAL] == 5.6
    @test costs[JADE.TimePoint(2019, 2)][:GAS] == 6.61
    @test costs[JADE.TimePoint(2019, 3)][:DIESEL] == 19.57
    @test costs[JADE.TimePoint(2019, 3)][:CO2] == 27.0
    @test fuels == Dict(:COAL => 0.0912, :DIESEL => 0.0, :GAS => 0.0528)
    @test length(costs) == 3
    return
end

function test_data_thermal_fuel_costs_non_contig()
    filename = _validation_file("thermal_fuel_costs_small_non_contig.csv")
    @test_throws(
        ErrorException("Weeks in $filename must be contiguous"),
        JADE.getfuelcosts(filename),
    )
    return
end

function test_data_reservoirs()
    reservoirs = JADE.initialisereservoirs(
        _validation_file("reservoirs_small.csv"),
        _validation_file("reservoir_limits_small.csv"),
    )
    @test length(reservoirs) == 2
    hawea = reservoirs[:LAKE_HAWEA]
    @test hawea.sp == 0.0
    @test hawea.index == 1
    @test hawea.initial == 332.074
    @test hawea.capacity isa JADE.TimeSeries{Float64}
    @test hawea.capacity.data == fill(1141.95, 52)
    @test hawea.contingent.data == [[JADE.ContingentTranche(0.0, 0.0)] for _ in 1:52]
    tekapo = reservoirs[:LAKE_TEKAPO]
    @test tekapo.initial == 397.0295
    @test tekapo.sp == 0.0
    @test tekapo.index == 2
    @test length(tekapo.capacity) == 52
    @test extrema(tekapo.capacity.data) == (514.1, 632.4)
    @test tekapo.contingent isa JADE.TimeSeries{Vector{JADE.ContingentTranche}}
    @test tekapo.contingent[1] ==
          [JADE.ContingentTranche(100.0, 500.0), JADE.ContingentTranche(91.0, 5000.0)]
    @test tekapo.contingent[30] ==
          [JADE.ContingentTranche(100.0, 0.0), JADE.ContingentTranche(91.0, 0.0)]
    return
end

function test_data_reservoirs_fail_duplicate()
    @test_throws(
        ErrorException("Reservoir LAKE_HAWEA given twice."),
        JADE.initialisereservoirs(
            _validation_file("reservoirs_fail_duplicate.csv"),
            _input_file("reservoir_limits.csv"),
        ),
    )
    return
end

function test_data_reservoirs_fail_nonconstant_min_level()
    @test_throws(
        ErrorException(
            "The maximum contingent storage is not constant for reservoir LAKE_TEKAPO.",
        ),
        JADE.initialisereservoirs(
            _validation_file("reservoirs_small.csv"),
            _validation_file("reservoir_limits_fail_nonconstant.csv"),
        ),
    )
    return
end

function test_data_reservoirs_fail_nonconstant_missing_penalty()
    @test_throws(
        ErrorException("LAKE_TEKAPO contingent storage missing MIN_2_PENALTY"),
        JADE.initialisereservoirs(
            _validation_file("reservoirs_small.csv"),
            _validation_file("reservoir_limits_fail_missing_penalty.csv"),
        ),
    )
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
    arcs = JADE.getnaturalarcs(_validation_file("hydro_arcs_small.csv"))
    @test length(arcs) == 2
    @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)] == JADE.NaturalArc(0.0, 100.0, -1.0, -1.0)
    @test arcs[(:LAKE_HAWEA, :LAKE_DUNSTAN)] == JADE.NaturalArc(10.0, Inf, -1.0, -1.0)
    return
end

function test_data_hydro_arcs_penalty()
    arcs = JADE.getnaturalarcs(_validation_file("hydro_arcs_small_penalty.csv"))
    @test length(arcs) == 2
    @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)] == JADE.NaturalArc(0.0, Inf, 100.0, -1.0)
    @test arcs[(:LAKE_HAWEA, :LAKE_DUNSTAN)] == JADE.NaturalArc(0.0, Inf, 25.0, 50.0)
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

function test_data_hydro_stations_small()
    stations, station_arcs =
        JADE.gethydros(_validation_file("hydro_stations_small.csv"), [:NI])
    @test length(stations) == 1
    @test length(station_arcs) == 1
    arc = (:LAKE_ARAPUNI, :LAKE_KARAPIRO)
    @test station_arcs[arc] == JADE.StationArc(Inf, -1.0, :ARAPUNI)
    @test stations[:ARAPUNI] == JADE.HydroStation(:NI, 196.7, 0.439847649, 0.0, arc)
    return
end

function test_data_hydro_stations_small_om_cost()
    stations, station_arcs =
        JADE.gethydros(_validation_file("hydro_stations_small_om_cost.csv"), [:NI])
    @test length(stations) == 1
    @test length(station_arcs) == 1
    arc = (:LAKE_ARAPUNI, :LAKE_KARAPIRO)
    @test station_arcs[arc] == JADE.StationArc(10, -1.0, :ARAPUNI)
    @test stations[:ARAPUNI] == JADE.HydroStation(:NI, 196.7, 0.439847649, 101.0, arc)
    return
end

function test_data_hydro_stations_small_penalty()
    stations, station_arcs =
        JADE.gethydros(_validation_file("hydro_stations_small_penalty.csv"), [:NI])
    @test length(stations) == 1
    @test length(station_arcs) == 1
    arc = (:LAKE_ARAPUNI, :LAKE_KARAPIRO)
    @test station_arcs[arc] == JADE.StationArc(Inf, 102.0, :ARAPUNI)
    @test stations[:ARAPUNI] == JADE.HydroStation(:NI, 196.7, 0.439847649, 101.0, arc)
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
