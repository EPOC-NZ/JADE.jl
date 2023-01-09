#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

"""
    JADEsddp(d::JADEData, optimizer = nothing)

This function builds a JADE model using SDDP.jl from the JADEData object `d`.
Returns an `SDDPModel` object.

### Required Arguments
`d` is a `JADEData` object.
`optimizer` is a `JuMP` optimizer.
"""
function JADEsddp(d::JADEData, optimizer = nothing)

    # Constants
    penalty_ub = d.rundata.penalty_ub
    penalty_lb = d.rundata.penalty_lb
    number_of_wks = d.rundata.number_of_wks
    nscenarios = d.rundata.nscenarios
    nmargins = length(d.terminal_eqns)
    scale_factor = d.rundata.scale_reservoirs
    scale_obj = d.rundata.scale_objective

    @assert nmargins > 0

    # Some more convenient names
    LOOPS = d.loops
    s = d.sets

    if optimizer == nothing
        error("No solver specified")
    elseif typeof(optimizer) <: Function
        d.parallel_optimizer = optimizer
        optimizer = d.parallel_optimizer()
    end

    #------------------------------------------------------------------------
    graph = SDDP.LinearGraph(number_of_wks)

    if d.rundata.steady_state
        SDDP.add_edge(graph, number_of_wks => 1, d.rundata.discount)
    end

    if d.rundata.weekly_discounting && d.rundata.discount != 0
        wk_disc = d.rundata.discount^(1 / 52)
        for (index, obj) in graph.nodes
            graph.nodes[index] = [(obj[1][1], wk_disc)]
        end
    end

    sddpm = SDDP.PolicyGraph(
        graph,
        sense = :Min,
        lower_bound = 0,
        optimizer = optimizer,
    ) do md, stage

        #-------------------------------------------------------------------------
        # Year and week for the current stage
        #-------------------------------------------------------------------------
        timenow = TimePoint(d.rundata.start_yr, d.rundata.start_wk) + stage - 1
        CONTINGENT = [
            r for r in s.RESERVOIRS if sum(
                d.reservoirs[r].contingent[timenow][j].level for
                j in 1:length(d.reservoirs[r].contingent[timenow])
            ) > 0.0
        ]

        dr_keys = [
            (n, bl) for n in keys(d.dr_tranches[timenow]), bl in s.BLOCKS if
            bl in keys(d.dr_tranches[timenow][n])
        ]

        en_keys = unique([
            (n, sector, lb) for n in keys(d.en_tranches[timenow]) for
            (sector, lb) in keys(d.en_tranches[timenow][n])
        ])

        #------------------------------------------------------------------------
        # State variable: water in reservoirs
        #------------------------------------------------------------------------
        JuMP.@variable(
            md,
            -sum(
                d.reservoirs[r].contingent[timenow][j].level / scale_factor for
                j in 1:length(d.reservoirs[r].contingent[timenow])
            ) / scale_factor <=
            reslevel[r in s.RESERVOIRS] <=
            d.reservoirs[r].capacity[timenow] / scale_factor,
            SDDP.State,
            initial_value = d.reservoirs[r].initial / scale_factor
        )

        #------------------------------------------------------------------------
        # Other variables
        #------------------------------------------------------------------------
        FLOWOVER = [a for a in s.NATURAL_ARCS if d.natural_arcs[a].maxflow != Inf]
        FLOWUNDER = [a for a in s.NATURAL_ARCS if d.natural_arcs[a].minflow != 0.0]
        SPILLOVER = [a for a in s.STATION_ARCS if d.station_arcs[a].maxflow != Inf]

        JuMP.@variables(
            md,
            begin
                # Dispatch of energy in MW from hydro stations
                hydro_disp[s.HYDROS, s.BLOCKS] >= 0
                # Amount of thermal energy used, in MW
                thermal_use[s.THERMALS, s.BLOCKS] >= 0
                # Transmission flows between nodes in MW
                transflow[s.TRANS_ARCS, s.BLOCKS]
                # Water flows in cumecs
                naturalflows[s.NATURAL_ARCS, s.BLOCKS] >= 0
                # Water going through hydro stations
                releases[s.STATION_ARCS, s.BLOCKS] >= 0
                # Water spilled from hydro stations
                spills[s.STATION_ARCS, s.BLOCKS] >= 0
                # Residual cost the final time period
                terminalcost >= 0
                # Lost load amount in MW
                0 <=
                lostload[
                    n in keys(d.dr_tranches[timenow]),
                    bl in keys(d.dr_tranches[timenow][n]),
                    k in keys(d.dr_tranches[timenow][n][bl]),
                ] <=
                d.dr_tranches[timenow][n][bl][k].q
                # Energy shedding
                0 <=
                energyshedding[
                    n in keys(d.en_tranches[timenow]),
                    s in keys(d.en_tranches[timenow][n]),
                    k in 1:length(d.en_tranches[timenow][n][s]),
                ] <=
                d.en_tranches[timenow][n][s][k].q * 1E3
                # Any flows over upper bound
                flowover[FLOWOVER, s.BLOCKS] >= 0
                # Any flows below lower bound
                flowunder[FLOWUNDER, s.BLOCKS] >= 0
                # Spill flows over upper bound (new)
                spillover[SPILLOVER, s.BLOCKS] >= 0
                # Contingent storage tranche
                contingent[
                    r in CONTINGENT,
                    1:length(d.reservoirs[r].contingent[timenow]),
                ] >= 0
                # To track inflow levels seen
                inflow[[s.CATCHMENTS_WITH_INFLOW; [:scenario]]]
            end
        )

        if d.rundata.losses != :none
            JuMP.@variables(
                md,
                begin
                    # Positive partof transflow
                    postransflow[s.TRANS_ARCS, s.BLOCKS] >= 0
                    # Negative part of transflow
                    negtransflow[s.TRANS_ARCS, s.BLOCKS] >= 0
                    # tranches of transmission flow for losses
                    postransflowtranche[
                        (i, j) in s.TRANS_ARCS,
                        1:length(d.transmission[(i, j)].Poslosses),
                        s.BLOCKS,
                    ] >= 0
                    negtransflowtranche[
                        (i, j) in s.TRANS_ARCS,
                        1:length(d.transmission[(i, j)].Neglosses),
                        s.BLOCKS,
                    ] >= 0
                    # losses in MW
                    losses[s.TRANS_ARCS, s.BLOCKS] >= 0
                    node_losses[s.NODES, s.BLOCKS] >= 0
                end
            )
        end

        #------------------------------------------------------------------------
        # Define handy expressions
        #------------------------------------------------------------------------

        JuMP.@expressions(
            md,
            begin
                # Number of hours in a week
                totHours, sum(d.durations[timenow][bl] for bl in s.BLOCKS)

                # Net transmission to any node: transmission to, minus transmission away
                transmission[n in s.NODES, bl in s.BLOCKS],
                sum(transflow[(i, j), bl] for (i, j) in s.TRANS_ARCS if j == n) -
                sum(transflow[(i, j), bl] for (i, j) in s.TRANS_ARCS if i == n)

                # Penalties for going over/under flow bounds. Note spMax is in MWh/m^3.
                flowpenalties,
                SECONDSPERHOUR * sum(
                    (
                        sum(
                            penalty_ub * d.spMax * flowover[a, bl] for
                            a in FLOWOVER if d.natural_arcs[a].ub_penalty == -1.0
                        ) +
                        sum(
                            penalty_ub * d.spMax * spillover[a, bl] for
                            a in SPILLOVER if d.station_arcs[a].penalty == -1.0
                        ) +
                        sum(
                            penalty_lb * d.spMax * flowunder[a, bl] for
                            a in FLOWUNDER if d.natural_arcs[a].lb_penalty == -1.0
                        ) +
                        (
                            sum(
                                d.natural_arcs[a].ub_penalty * flowover[a, bl] for
                                a in FLOWOVER if d.natural_arcs[a].ub_penalty > 0
                            ) +
                            sum(
                                d.station_arcs[a].penalty * spillover[a, bl] for
                                a in SPILLOVER if d.station_arcs[a].penalty > 0
                            ) +
                            sum(
                                d.natural_arcs[a].lb_penalty * flowunder[a, bl] for
                                a in FLOWUNDER if d.natural_arcs[a].lb_penalty > 0
                            )
                        ) / 1000
                    ) * d.durations[timenow][bl] for bl in s.BLOCKS
                )

                # Flow in minus flow out to any node
                netflow[n in s.CATCHMENTS, bl in s.BLOCKS],
                sum(naturalflows[(i, j), bl] for (i, j) in s.NATURAL_ARCS if j == n) -
                sum(naturalflows[(i, j), bl] for (i, j) in s.NATURAL_ARCS if i == n) +
                sum(releases[(i, j), bl] for (i, j) in s.STATION_ARCS if j == n) -
                sum(releases[(i, j), bl] for (i, j) in s.STATION_ARCS if i == n) +
                sum(spills[(i, j), bl] for (i, j) in s.STATION_ARCS if j == n) -
                sum(spills[(i, j), bl] for (i, j) in s.STATION_ARCS if i == n)
            end
        )

        # Total supply of electricity at any node and block
        if d.rundata.losses != :none
            JuMP.@expression(
                md,
                supply[n in s.NODES, bl in s.BLOCKS],
                d.durations[timenow][bl] * (
                    transmission[n, bl] - node_losses[n, bl] +
                    sum(thermal_use[m, bl] for m in d.nodehas[n].thermal) +
                    sum(hydro_disp[m, bl] for m in d.nodehas[n].hydro)
                )
            )
        else
            JuMP.@expression(
                md,
                supply[n in s.NODES, bl in s.BLOCKS],
                d.durations[timenow][bl] * (
                    transmission[n, bl] +
                    sum(thermal_use[m, bl] for m in d.nodehas[n].thermal) +
                    sum(hydro_disp[m, bl] for m in d.nodehas[n].hydro)
                )
            )
        end

        #------------------------------------------------------------------------
        # Define constraints
        #------------------------------------------------------------------------
        JuMP.@constraints(
            md,
            begin
                # Lower and upper bounds on flows
                natOver[a in FLOWOVER, bl in s.BLOCKS],
                flowover[a, bl] >= naturalflows[a, bl] - d.natural_arcs[a].maxflow
                natUnder[a in FLOWUNDER, bl in s.BLOCKS],
                flowunder[a, bl] >= d.natural_arcs[a].minflow - naturalflows[a, bl]
                spillOver[a in SPILLOVER, bl in s.BLOCKS],
                spillover[a, bl] >= spills[a, bl] - d.station_arcs[a].maxflow

                # Capacity constraints

                # Hydro plant capacities
                useHydro[m in s.HYDROS, bl in s.BLOCKS],
                hydro_disp[m, bl] <=
                d.hydro_stations[m].capacity - sum(
                    d.outage[timenow][(mm, bb)] for
                    (mm, bb) in keys(d.outage[timenow]) if (mm, bb) == (m, bl)
                )

                # Thermal plant capacities
                useThermal[m in s.THERMALS, bl in s.BLOCKS],
                thermal_use[m, bl] <=
                d.thermal_stations[m].capacity - sum(
                    d.outage[timenow][(mm, bb)] for
                    (mm, bb) in keys(d.outage[timenow]) if (mm, bb) == (m, bl)
                )

                # Transmission line capacities
                transUpper[(n, m) in s.TRANS_ARCS, bl in s.BLOCKS],
                transflow[(n, m), bl] <=
                d.transmission[(n, m)].poscapacity -
                d.transmission[(n, m)].posoutage[timenow][bl]
                transLower[(n, m) in s.TRANS_ARCS, bl in s.BLOCKS],
                -d.transmission[(n, m)].negcapacity +
                d.transmission[(n, m)].negoutage[timenow][bl] <= transflow[(n, m), bl]

                # Set thermal station capacities to zero if the station has not yet been
                #    commissioned, or is decommissioned
                notCommissioned[
                    t in s.THERMALS,
                    bl in s.BLOCKS;
                    timenow < d.thermal_stations[t].commission,
                ],
                thermal_use[t, bl] <= 0
                decommissioned[
                    t in s.THERMALS,
                    bl in s.BLOCKS;
                    d.thermal_stations[t].decommission != TimePoint(0, 0) &&
                        timenow > d.thermal_stations[t].decommission,
                ],
                thermal_use[t, bl] <= 0

                # All flow through hydro-stations is used in hydro generation

                defineDispatch[m in s.HYDROS, bl in s.BLOCKS],
                hydro_disp[m, bl] ==
                d.hydro_stations[m].sp * releases[d.hydro_stations[m].arc, bl]

                # Define shedding: the load shed over all sectors has to equal the shortage. In MWh.

                defineShedding[n in s.NODES, bl in s.BLOCKS],
                sum(lostload[n, bl, k] for k in keys(d.dr_tranches[timenow][n][bl])) *
                d.durations[timenow][bl] >=
                d.demand[timenow][(n, bl)] - d.fixed[timenow][(n, bl)] - supply[n, bl]

                energyShedding[(n, sector, loadblocks) in en_keys],
                sum(
                    lostload[n, bl, (s, name)] * d.durations[timenow][bl] for
                    bl in loadblocks,
                    (s, name) in keys(d.dr_tranches[timenow][n][bl]) if s == sector
                ) <= sum(
                    energyshedding[n, (sector, loadblocks), k] for
                    k in 1:length(d.en_tranches[timenow][n][(sector, loadblocks)])
                )
            end
        )

        if length(CONTINGENT) != 0
            JuMP.@constraints(
                md,
                begin
                    contingentstorage[r in CONTINGENT],
                    reslevel[r].out >=
                    -sum(
                        contingent[r, j] / scale_factor for
                        j in 1:length(d.reservoirs[r].contingent[timenow])
                    )

                    maxcontingenttranche[
                        r in CONTINGENT,
                        j in 1:(length(d.reservoirs[r].contingent[timenow])-1),
                    ],
                    contingent[r, j] <= d.reservoirs[r].contingent[timenow][j].level
                end
            )
        end

        ###################
        #ABP losses code
        if d.rundata.losses != :none
            JuMP.@constraints(
                md,
                begin
                    # Define flow tranches for piecewise linear losses.
                    definePosFlowTranches[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    postransflow[(i, j), bl] == sum(
                        postransflowtranche[(i, j), k, bl] for
                        k in 1:length(d.transmission[(i, j)].Poslosses)
                    )

                    defineNegFlowTranches[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    negtransflow[(i, j), bl] == sum(
                        negtransflowtranche[(i, j), k, bl] for
                        k in 1:length(d.transmission[(i, j)].Neglosses)
                    )

                    definePosLossTranche1[
                        (i, j) in s.TRANS_ARCS,
                        bl in s.BLOCKS,
                        k in 1:length(d.transmission[(i, j)].Poslosses),
                    ],
                    postransflowtranche[(i, j), k, bl] <=
                    d.transmission[(i, j)].Poslosses[k][1]

                    defineNegLossTranche1[
                        (i, j) in s.TRANS_ARCS,
                        bl in s.BLOCKS,
                        k in 1:length(d.transmission[(i, j)].Neglosses),
                    ],
                    negtransflowtranche[(i, j), k, bl] <=
                    d.transmission[(i, j)].Neglosses[k][1]

                    defineArcLosses[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    losses[(i, j), bl] >=
                    sum(
                        postransflowtranche[(i, j), k, bl] *
                        d.transmission[(i, j)].Poslosses[k][2] for
                        k in 1:length(d.transmission[(i, j)].Poslosses)
                    ) + sum(
                        negtransflowtranche[(i, j), k, bl] *
                        d.transmission[(i, j)].Neglosses[k][2] for
                        k in 1:length(d.transmission[(i, j)].Neglosses)
                    )

                    # Define positive and Negative part of transflow
                    definePosFlow[(i, j) in s.TRANS_ARCS, bl in s.BLOCKS],
                    transflow[(i, j), bl] ==
                    postransflow[(i, j), bl] - negtransflow[(i, j), bl]
                end
            )

            # Half line losses are to be added to load at each end of the line
            JuMP.@constraints(
                md,
                begin
                    defineNodalLosses[n in s.NODES, bl in s.BLOCKS],
                    node_losses[n, bl] ==
                    0.5 * sum(losses[(i, j), bl] for (i, j) in s.TRANS_ARCS if j == n) +
                    0.5 * sum(losses[(i, j), bl] for (i, j) in s.TRANS_ARCS if i == n)
                end
            )
        end

        # Power times reactance of arcs in a loop adds to zero
        # ABP disable
        #       loopflow[l in LOOPS, bl in s.BLOCKS],
        #           sum(d.transmission[a].reactance * transflow[a, bl] for a in l) == 0
        #     end)

        #ABP losses code ends
        #####################

        #------------------------------------------------------------------------
        # Flow-conservation-related calculations
        #------------------------------------------------------------------------
        # Get an inflow for every location with inflow data
        # JuMP.@variable(md, ω[s.CATCHMENTS_WITH_INFLOW])
        # for c in s.CATCHMENTS_WITH_INFLOW
        #     JuMP.@constraint(md, inflow[c] == ω[c])
        # end

        inflow_uncertainty = Array{Dict{Symbol,Float64},1}()
        for scenario in 1:length(d.inflow_mat[1][s.CATCHMENTS_WITH_INFLOW[1]])
            s_inflows = Dict{Symbol,Float64}()
            for c in s.CATCHMENTS_WITH_INFLOW
                s_inflows[c] = d.inflow_mat[timenow.week][c][scenario]
            end
            if stage == 1 && d.rundata.first_week_known
                s_inflows[:scenario] = 0
            else
                s_inflows[:scenario] = scenario
            end
            push!(inflow_uncertainty, s_inflows)
        end
        SDDP.parameterize(md, inflow_uncertainty) do ϕ
            for (c, value) in ϕ
                JuMP.fix(inflow[c], value)
            end
        end

        JuMP.@constraints(
            md,
            begin
                # Conservation for reservoirs
                rbalance[r in s.RESERVOIRS],
                (reslevel[r].out - reslevel[r].in) * 1E3 * scale_factor ==
                SECONDSPERHOUR / 1E3 * (
                    sum(d.durations[timenow][bl] * (netflow[r, bl]) for bl in s.BLOCKS) + totHours * inflow[r]
                )

                # Conservation for junction points with inflow
                jbalance[c in s.CATCHMENTS_WITH_INFLOW, bl in s.BLOCKS; c in s.JUNCTIONS],
                netflow[c, bl] + inflow[c] == 0

                # Flow conservation for junctions without an inflow
                conserveFlow[c in s.JUNCTIONS_WITHOUT_INFLOW, bl in s.BLOCKS],
                netflow[c, bl] == 0
            end
        )

        for dr in d.rundata.decision_rules
            if timenow.week ∉ dr.weeks
                continue
            end
            LHS = 0.0
            if dr.flowtype == :generation
                LHS =
                    d.hydro_stations[dr.station].sp * sum(
                        releases[d.hydro_stations[dr.station].arc, bl] *
                        d.durations[timenow][bl] for bl in s.BLOCKS
                    )
            elseif dr.flowtype == :spill
                LHS =
                    d.hydro_stations[dr.station].sp * sum(
                        spills[d.hydro_stations[dr.station].arc, bl] *
                        d.durations[timenow][bl] for bl in s.BLOCKS
                    )
            elseif dr.flowtype == :combined
                LHS =
                    d.hydro_stations[dr.station].sp * sum(
                        (
                            releases[d.hydro_stations[dr.station].arc, bl] +
                            spills[d.hydro_stations[dr.station].arc, bl]
                        ) * d.durations[timenow][bl] for bl in s.BLOCKS
                    )
            else
                error("Invalid flow type: " * string(dr.flowtype))
            end
            RHS =
                dr.intercept +
                dr.slope * (
                    reslevel[dr.reservoir].in * scale_factor +
                    SECONDSPERHOUR * totHours * inflow[dr.reservoir] / 1E6
                )

            if dr.boundtype == :upper
                JuMP.@constraint(md, LHS <= RHS)
            elseif dr.boundtype == :lower
                JuMP.@constraint(md, LHS >= RHS)
            elseif dr.boundtype == :equality
                JuMP.@constraint(md, LHS == RHS)
            else
                error("Invalid bound type: " * string(dr.boundtype))
            end
        end

        #------------------------------------------------------------------------
        # Objective-related calculations
        #------------------------------------------------------------------------
        JuMP.@expression(
            md,
            lostloadcosts,
            sum(
                sum(
                    lostload[n, bl, k] *
                    d.dr_tranches[timenow][n][bl][k].p *
                    d.durations[timenow][bl] for
                    k in keys(d.dr_tranches[timenow][n][bl])
                ) for (n, bl) in dr_keys
            ) + sum(
                sum(
                    energyshedding[n, (sector, loadblocks), k] *
                    d.en_tranches[timenow][n][(sector, loadblocks)][k].p for
                    k in 1:length(d.en_tranches[timenow][n][(sector, loadblocks)])
                ) for (n, sector, loadblocks) in en_keys
            )
        )

        JuMP.@expression(
            md,
            contingent_storage_cost,
            sum(
                contingent[r, j] / scale_factor *
                d.reservoirs[r].contingent[timenow][j].penalty for r in CONTINGENT,
                j in 1:length(d.reservoirs[r].contingent[timenow])
            )
        )

        JuMP.@expression(
            md,
            carbon_emissions[t in s.THERMALS, bl in s.BLOCKS],
            d.carbon_content[d.thermal_stations[t].fuel] *
            d.thermal_stations[t].heatrate *
            thermal_use[t, bl] *
            d.durations[timenow][bl]
        )

        JuMP.@expression(
            md,
            immediate_cost,
            sum(
                (station.omcost + d.fuel_costs[timenow][station.fuel] * station.heatrate) *
                thermal_use[name, bl] *
                d.durations[timenow][bl] +
                carbon_emissions[name, bl] * d.fuel_costs[timenow][:CO2] for
                (name, station) in d.thermal_stations, bl in s.BLOCKS
            ) +
            sum(
                station.omcost * hydro_disp[name, bl] * d.durations[timenow][bl] for
                (name, station) in d.hydro_stations, bl in s.BLOCKS
            ) +
            flowpenalties +
            lostloadcosts +
            contingent_storage_cost
        )

        if stage < number_of_wks || !d.rundata.use_terminal_mwvs
            # Stage cost function not including terminal water value
            SDDP.@stageobjective(md, immediate_cost / scale_obj)
        else
            # Convert stored water in Mm³ to MWh
            JuMP.@expression(
                md,
                storedenergy,
                1E6 *
                scale_factor *
                sum(d.reservoirs[r].sp * reslevel[r].out for r in s.RESERVOIRS)
            )

            for cut in d.terminal_eqns # -terminalcost for value
                JuMP.@constraint(
                    md,
                    (-terminalcost) <=
                    (cut.intercept + cut.coefficient * storedenergy) / scale_obj
                )
            end
            # Cost function includes terminal values added
            SDDP.@stageobjective(md, immediate_cost / scale_obj + terminalcost)
        end
    end

    return sddpm
end
