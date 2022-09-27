function inputfile(s)
    return joinpath(JADE.@JADE_DIR, "Input", "test1", s)
end

TimePoint = JADE.TimePoint

@testset "Data Input Validation" begin
    @testset "Thermal Stations" begin
        stations =
            JADE.getthermalstations(inputfile("thermal_stations.csv"), [:NI, :HAY, :SI])
        @test length(stations) == 13
        stratford = stations[:STRATFORD_PEAKERS]
        @test stratford.fuel == :GAS
        @test stratford.heatrate == 9.5
    end

    @testset "Costs" begin
        costs, fuels = JADE.getfuelcosts(inputfile("thermal_fuel_costs.csv"))
        @test length(fuels) == 3
        @test costs[TimePoint(2005, 10)][:COAL] == 4.0
        @test costs[TimePoint(2009, 52)][:DIESEL] == 24.33024581
        @test length(costs) == 260
    end

    @testset "Reservoirs" begin
        reservoirs = JADE.initialisereservoirs(
            inputfile("reservoirs.csv"),
            inputfile("reservoir_limits.csv"),
        )
        @test length(reservoirs) == 7
        benmore = reservoirs[:LAKE_BENMORE]
        @test benmore.initial == 322.00032233399
    end

    @testset "Transmission Lines" begin
        lineoutage, arcs = JADE.gettimeseries(inputfile("transmission_outages.csv"))
        lineoutage = JADE.duplicate_data(lineoutage, [:PEAK, :SHOULDER, :OFFPEAK])

        arcs, losses = JADE.gettransarcs(
            inputfile("transmission.csv"),
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
    end

    @testset "Natural Arcs" begin
        arcs = JADE.getnaturalarcs(inputfile("hydro_arcs.csv"))
        @test length(arcs) == 17
        @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)].minflow == 0.0
        @test arcs[(:LAKE_WANAKA, :LAKE_DUNSTAN)].maxflow == Inf
    end

    @testset "Hydro Stations" begin
        stations, station_arcs =
            JADE.gethydros(inputfile("hydro_stations.csv"), [:NI, :HAY, :SI])
        @test length(stations) == 26
        tokaanu = stations[:TOKAANU]
        @test tokaanu.node == :NI
        @test tokaanu.capacity == 240.0
    end

    @testset "Outages" begin
        outages, stations = JADE.gettimeseries(inputfile("station_outages.csv"))
        @test length(outages) == 260
    end

    @testset "Terminal" begin
        eqn = JADE.getterminalvalue(inputfile("terminal_water_value.csv"))
        @test length(eqn) == 5
        @test eqn[1].coefficient == 137.218398630355
    end

    @testset "Demand / Durations" begin
        duration, blocks = JADE.gettimeseries(inputfile("hours_per_block.csv"))
        demand, nodes = JADE.getdemand(inputfile("demand.csv"), duration)
        @test length(nodes) == 3
        @test setdiff(nodes, [:NI, :SI, :HAY]) == []
        @test length(demand) == 260
        @test length(keys(demand[TimePoint(2005, 1)])) == 9
        @test isapprox(demand[TimePoint(2005, 1)][(:SI, :OFFPEAK)], 60764.5, atol = 0.1)
        @test length(blocks) == 3
        @test setdiff(blocks, [:OFFPEAK, :PEAK, :SHOULDER]) == []
        @test length(duration) == 260
        @test length(keys(duration[TimePoint(2005, 1)])) == 3
        @test duration[TimePoint(2005, 1)][:OFFPEAK] == 51
    end

    # @testset "Durations" begin end

    # @testset "Shedding" begin
    #     lost_load_cost, sectors = JADE.getshedding(inputfile("lost_load.csv"))
    #     @test length(sectors) == 3
    #     @test setdiff(sectors, [:INDUSTRIAL, :COMMERCIAL, :RESIDENTIAL]) == []
    #     @test length(lost_load_cost[(:NI, :INDUSTRIAL)]) == 3
    # end
end
