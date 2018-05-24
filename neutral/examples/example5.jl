# Configuration for running simple neutral copy error simulation
#=
Recommended command line to run:
> julia -p 4 -L NeutralEvolution.jl nrun.jl examples/example1
> julia -p 4 -L NeutralEvolution.jl nrun.jl examples/example1 1  # include 1 to set random number seed
=#
global simtype
@everywhere simtype = 3
export simtype
@everywhere const N = 100        # population size
const mutstddev = 0.04
const ngens = 351
const initial_value = 1.0
const num_trials=1000
const record_interval = 350
const use_population = true  
#const use_population = false
const fit_slope=1.0
const log_error = false
const wright_fisher_copy = false
const neutral = true
