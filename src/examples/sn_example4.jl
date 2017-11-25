# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
>  julia neutral.jl examples/sn_example4 1   # include 1 to set random number seed
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 1        # population size
const mutstddev = 0.04
const ngens = 10
const initial_value = 1.0
const num_trials=1
const record_interval = 1
const use_population = true  # try this both with and without
#const use_population = false
const log_error = true
