# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example2
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 10        # population size
const mutstddev = 0.04
const ngens = 101
const initial_value = 1.0
const num_trials=8
const record_interval = 10
const use_population=true
const save_populations=true
