# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example1
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 10        # population size
const mutstddev = 0.04
const ngens = 1001
const initial_value = 1.0
const num_trials=1000
const record_interval = 200
const log_error = false
const wright_fisher_copy = false
const conformist_probability=0.3
