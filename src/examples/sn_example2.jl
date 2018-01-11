# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example2
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 5        # population size
const mutstddev = 0.04
const ngens = 6
const initial_value = 1.0
const num_trials=3
const record_interval = 2
const use_population=true
const save_populations=true
const conformist_probability=1.0
