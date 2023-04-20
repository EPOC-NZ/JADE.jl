#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

using Test

for file in filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(@__DIR__))
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
