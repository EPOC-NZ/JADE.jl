#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

using JuliaFormatter

function format_code()
    i = 0
    while i < 10
        if format(dirname(@__DIR__))
            if i == 0
                return "Formatting unchanged"
            else
                return "Formatting corrected"
            end
        end
        i += 1
    end
    return "Formatting failed"
end

@info("Formatting...")
@info(format_code())
