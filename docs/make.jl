#  This file is part of JADE source code.
#  Copyright © 2016-2022 Electric Power Optimization Centre, University of Auckland.
#
#  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
#  If a copy of the MPL was not distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

using Documenter, JADE

makedocs(
    sitename = "JADE: a Julia DOASA Environment",
    modules = [JADE],
    clean = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["JADE" => "index.md", "API Reference" => "api.md"],
)

deploydocs(repo = "github.com/EPOC-NZ/JADE.git", devurl = "docs")
