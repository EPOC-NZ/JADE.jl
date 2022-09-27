# User interface functions
These functions are to be utilised by JADE users to set up and run JADE models.

## Functions to set up JADE model / data

### Using run files
```@docs
JADE.define_JADE_model(::String)
JADE.define_JADE_solve_options(::String)
JADE.define_JADE_simulation(::String)
```

### Using Julia scripts
```@docs
JADE.define_JADE_model()
JADE.define_JADE_solve_options()
JADE.define_JADE_simulation()
JADE.create_JADE_model
JADE.DecisionRule
```

## Functions for training and simulation
```@docs
JADE.optimize_policy!
JADE.simulate
```

## Visualisation functions
```@docs
JADE.plot_storage
JADE.plot_prices
```

# Internal functions
These are internal functions within JADE that should only be used or modified
by persons who have a deep understanding of the underlying SDDP model.

## Utility functions
```@docs
JADE.gettimeseries
JADE.calculatetime
JADE.cleandata
JADE.between
```

## Data input functions
```@docs
JADE.parsefile
JADE.getitems
JADE.parse_run_options
JADE.gethydros
JADE.getinflows_for_historical
JADE.getnaturalarcs
JADE.getterminalvalue
JADE.getthermalstations
JADE.initialisereservoirs
JADE.getdemand
JADE.checkoutages
JADE.get_file_directory
JADE.check_settings_compatibility
```

## Data processing functions   
```@docs
JADE.out_neighbors
JADE.hasdownstream
JADE.adjustinflows
```

## Data output functions
```@docs
JADE.diatofile
JADE.write_sim_results
JADE.write_training_results
```

## Structs
```@docs
JADE.RunData
JADE.JADESolveOptions
JADE.JADESimulation
JADE.JADEData
```

## Modified SDDP.jl functionality
```@docs
JADE.read_finalcuts_from_file
JADE.write_cuts_to_file
JADE.WrapHistorical
JADE.InSampleMonteCarlo2
```
