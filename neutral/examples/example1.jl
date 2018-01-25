# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example1
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 3        # population size
const mutstddev = 0.04
const ngens = 2
const initial_value = 1.0
const num_trials=1
const record_interval = 1
const log_error = false
const wright_fisher_copy = false
const conformist_probability=0.0
