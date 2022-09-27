#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

module JADE

using DataFrames, JuMP, SDDP, Random, DelimitedFiles, JSON

const SECONDSPERHOUR = 3600
const WEEKSPERYEAR = 52

macro JADE_DIR()
    ex = quote
        get(ENV, "JADE_DIR", "")::String
    end
    return esc(ex)
end

JSON.lower(t::SDDP.Expectation) = (1.0, 1.0)
function JSON.lower(t::T where {T<:SDDP.ConvexCombination})
    return (t.measures[1][1], t.measures[2][2].β)
end

include("structs.jl")
include("data.jl")                # data file of JADE
include("utilities.jl")
include("model.jl")               # model file of JADE
include("results.jl")             # functions for writing results to file
include("solve.jl")
include("simulate.jl")
include("visualise.jl")
include("sddp_modifications.jl")

export JADEdata,
    JADEmodel,
    optimize_policy!,
    define_JADE_model,
    define_JADE_solve_options,
    create_JADE_model,
    define_JADE_simulation,
    simulate,
    load_model_parameters,
    load_simulation_parameters
end
