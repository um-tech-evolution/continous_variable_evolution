# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia run.jl examples/example1
or
>  julia -p 4 run.jl examples/example1
=#
export simtype
simtype = 1 
const num_trials=1
const N_list = [6]        # Meta-population size list
#const N = 8        # Meta-population size
#const mutation_stddev = 0.05
const mutation_stddev_list = [0.14]
#const num_attributes = 2        # number attributes for quantitative representation
#const num_attributes_list = [2]        # number attributes for quantitative representation
const num_attributes_list = [1]        # number attributes for quantitative representation
#const ngens = 3       # Generations after burn-in
const ngens = 4       # Generations after burn-in
const num_subpops = 1
#const burn_in= 3    # if integer, generations of burn in.  If float, int_burn_in=burn_in*N+50
const burn_in= 1    # if integer, generations of burn in.  If float, int_burn_in=burn_in*N+50
const ideal = 1.0
const fit_slope = 1.0
const neutral = false
const w = 0.1
