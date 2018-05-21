# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
> julia -p 4 -L NeutralEvolution.jl nrun.jl examples/example1
=#
export simtype
@everywhere simtype = 3
@everywhere const N = 1        # population size
const mutstddev = 0.04
const ngens = 1000
#const initial_value = exp(5.0)
const initial_value = 1.0
const num_trials=100000
const record_interval = 50
const use_population = true
const log_error = true
const neutral = true
