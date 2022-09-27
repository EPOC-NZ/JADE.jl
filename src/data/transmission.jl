#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

#---------------------------------------------------
# Transmission information
#---------------------------------------------------

##############################################
#ABP transmission model changes here
##############################################
struct TransArc
    poscapacity::Float64
    negcapacity::Float64
    Poslosses::Vector{Tuple{Float64,Float64}}#Array{Float64,2}
    Neglosses::Vector{Tuple{Float64,Float64}}#Array{Float64,2}
    Presistance::Float64
    Nresistance::Float64
    posoutage::TimeSeries{Dict{Symbol,Float64}}
    negoutage::TimeSeries{Dict{Symbol,Float64}}
    # reactance::Float64
end

function gettransarcs(
    file::String,
    lineoutage::TimeSeries{Dict{Tuple{Symbol,Symbol},Float64}},
    losses::Symbol,
    loadblocks::Vector{Symbol},
)
    transmission = Dict{NTuple{2,Symbol},TransArc}()
    presistance = 0.0
    nresistance = 0.0
    vSPDlossArray = [
        0.14495 0.1087125
        0.32247 0.46742
        0.5 0.82247
        0.67753 1.17753
        0.85505 1.53258
        1.0 1.8912875
    ]

    poutage = TimeSeries{Float64}
    noutage = TimeSeries{Float64}

    first = true
    parsefile(file, true) do items
        if first
            if losses == :default
                @assert length(items) >= 3
                if length(items) == 3
                    losses = :none
                elseif length(items) == 4
                    losses = :quadratic
                else
                    losses = :piecewise
                end
            end
            @info("Losses mode: " * string(losses))
            if losses == :none
                @assert length(items) >= 3
            elseif losses == :piecewise
                if lowercase(items[1]) == "from_node"
                    @assert length(items) > 4
                    @assert lowercase(items[4]) == "loss_tranche_1"
                end
            elseif losses == :resistive
                if lowercase(items[1]) == "from_node"
                    @assert length(items) == 4
                    @assert lowercase(items[4]) == "resistance"
                end
            end
            first = false
        elseif losses ∈ [:piecewise, :none]
            # PIECEWISE LINEAR MODEL
            line_name = items[1] * "_TO_" * items[2]
            #			@info(" Line name = ",line_name)
            od_pair = (str2sym(items[1]), str2sym(items[2]))
            if haskey(transmission, od_pair)
                @warn("Transmission line $(od_pair) supplied before, so skipped.")
            else
                lossArray = Tuple{Float64,Float64}[]
                previous = 0.0
                i = 4

                while i < length(items) && losses == :piecewise
                    temp = parse(Float64, items[i]) - previous
                    if temp < 0.0
                        error(
                            "Loss tranches on transmission line $line_name must be increasing",
                        )
                    end
                    push!(lossArray, (temp, parse(Float64, items[i+1])))
                    previous = temp + previous
                    i = i + 2
                end
                if losses == :piecewise && i == 4
                    error(
                        "Line $line_name specified without any losses.\nIf there are no losses, specify a loss tranche with the line capacity and a 0 loss fraction.",
                    )
                end

                capacity = parse(Float64, items[3])

                if losses == :piecewise &&
                   sum(lossArray[j][1] for j in 1:length(lossArray)) < capacity
                    @warn(
                        "Loss tranches for line $line_name are less than capacity of line."
                    )
                end

                temp = Dict{Symbol,Float64}[]
                for j in 1:length(lineoutage)
                    temp2 = Dict{Symbol,Float64}()
                    for b in loadblocks
                        if haskey(lineoutage[j], (str2sym(line_name), b))
                            temp2[b] = lineoutage[j][(str2sym(line_name), b)]
                            delete!(lineoutage[j], (str2sym(line_name), b))
                            if temp2[b] > capacity
                                error(
                                    "An outage specified for line " *
                                    line_name *
                                    " is larger than the capacity of the line.",
                                )
                            end
                        else
                            temp2[b] = 0.0
                        end
                    end
                    push!(temp, temp2)
                end

                outage = TimeSeries{Dict{Symbol,Float64}}(lineoutage.startpoint, temp)

                if haskey(transmission, (od_pair[2], od_pair[1]))
                    @info(
                        "Transmission line $(od_pair) given for both directions, setting additional properties."
                    )

                    Pcapacity = transmission[(od_pair[2], od_pair[1])].poscapacity
                    poutage = transmission[(od_pair[2], od_pair[1])].posoutage
                    plossArray = transmission[(od_pair[2], od_pair[1])].Poslosses

                    transmission[(od_pair[2], od_pair[1])] = TransArc(
                        Pcapacity,
                        capacity,
                        plossArray,
                        lossArray,
                        0.0,
                        0.0,
                        poutage,
                        outage,
                    )
                else
                    transmission[od_pair] = TransArc(
                        capacity,
                        capacity,
                        lossArray,
                        lossArray,
                        0.0,
                        0.0,
                        outage,
                        outage,
                    )
                end
            end
        elseif losses == :resistive
            #  RESISTANCE MODEL
            # losses are 0.01*resistance*flow^2
            line_name = items[1] * "_TO_" * items[2]
            #			@info(" Line name = ",line_name)
            od_pair = (str2sym(items[1]), str2sym(items[2]))
            if haskey(transmission, od_pair)
                @warn("Transmission line $(od_pair) supplied before, so skipped.")
            else
                if haskey(transmission, (od_pair[2], od_pair[1]))
                    @info(
                        "Transmission line $(od_pair) given for both directions, setting additional properties."
                    )
                    plossArray = transmission[(od_pair[2], od_pair[1])].Poslosses
                    nlossArray = Tuple{Float64,Float64}[]
                    Ncapacity = parse(Float64, items[3])
                    Pcapacity = transmission[(od_pair[2], od_pair[1])].poscapacity
                    presistance = transmission[(od_pair[2], od_pair[1])].Presistance
                    nresistance = parse(Float64, items[4])
                    @info("presistance = ", presistance)
                    @info("nresistance = ", nresistance)
                    previous = 0.0
                    for k in 1:6
                        push!(
                            nlossArray,
                            (
                                vSPDlossArray[k, 1] * Ncapacity - previous,
                                vSPDlossArray[k, 2] * Ncapacity * 0.01 * nresistance,
                            ),
                        )
                        previous = vSPDlossArray[k, 1] * Ncapacity
                    end

                    temp = Float64[]
                    for j in 1:length(lineoutage)
                        push!(temp, lineoutage[j][str2sym(line_name)])
                    end
                    noutage = TimeSeries{Float64}(lineoutage.startpoint, temp)

                    transmission[(od_pair[2], od_pair[1])] = TransArc(
                        Pcapacity,
                        Ncapacity,
                        plossArray,
                        nlossArray,
                        presistance,
                        nresistance,
                        poutage,
                        noutage,
                    )
                else
                    #   This is the first time the node pair encountered so fill up the loss array and make poscapacity=negcapacity
                    plossArray = Tuple{Float64,Float64}[]
                    presistance = parse(Float64, items[4])
                    nresistance = presistance
                    Pcapacity = parse(Float64, items[3])
                    previous = 0.0

                    for k in 1:6
                        push!(
                            plossArray,
                            (
                                vSPDlossArray[k, 1] * Pcapacity - previous,
                                vSPDlossArray[k, 2] * Pcapacity * 0.01 * presistance,
                            ),
                        )
                        previous = vSPDlossArray[k, 1] * Pcapacity
                    end

                    temp = Float64[]
                    for j in 1:length(lineoutage)
                        push!(temp, lineoutage[j][str2sym(line_name)])
                    end
                    poutage = TimeSeries{Float64}(lineoutage.startpoint, temp)
                    transmission[od_pair] = TransArc(
                        Pcapacity,
                        Pcapacity,
                        plossArray,
                        plossArray,
                        presistance,
                        nresistance,
                        poutage,
                        poutage,
                    )
                end
            end
        end
    end

    for j in 1:length(lineoutage)
        if length(keys(lineoutage[j])) != 0
            error(
                "Invalid line name(s) specified in 'line_outages.csv' or 'transmission_outages.csv'.\n Only lines in 'transmission.csv' can have outages specified.",
            )
        end
    end

    return transmission, losses
end

##############################################
#ABP transmission model ends here
##############################################
