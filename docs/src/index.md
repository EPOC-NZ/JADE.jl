```@meta
CurrentModule = JADE
DocTestSetup = quote
    using JADE
end
```
# JADE: a Julia DOASA Environment
JADE is a sophisticated hydro-thermal (and demand-response) scheduling tool calibrated to the New
Zealand electricity system. JADE leverages the [SDDP.jl](https://github.com/odow/SDDP.jl) package
in order to find optimal waters values and release policies.

## Requirements

JADE requires [Julia](https://julialang.org/downloads/) 1.6+, JuMP 1.0+, SDDP.jl and appropriate optimizer(s).

## Installation

JADE is installed by the `pkg` interface provided by Julia. In the Julia REPL,
simply enter the following command.

    ] add "https://github.com/EPOC-NZ/JADE.git"

Julia will install all dependencies automatically, ensuring compatibility. If you encounter an issue
during installation, it may be due to another package with dependencies that conflict with JADE's.

## Running JADE

See the [documentation](https://github.com/EPOC-NZ/JADE/raw/main/docs/JADE%20documentation.pdf)
for specifics about the JADE inputs and how to run JADE.

The Electricity Authority hosts a repository of JADE input files. These can be accessed from the [EMI website](https://www.emi.ea.govt.nz/Wholesale/Tools/Jade).
Any input files in the JADE Github repository are for testing purposes only.

## Bugs

Please raise an [issue](https://github.com/EPOC-NZ/JADE/issues) if you experience an error while using JADE.
