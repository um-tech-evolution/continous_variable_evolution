# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example1
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 10        # population size
const mutstddev = 0.04
const ngens = 1000
const initial_value = 1.0
const num_trials=1000
const record_interval = 50
