# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run.jl examples/example1
or
>  julia -p 4 -L ContVarEvolution.jl run.jl examples/example1
=#
export simtype
@everywhere simtype = 2    
const num_trials=1
@everywhere const N_list = [3]        # Meta-population size list
#@everywhere const N = 8        # Meta-population size
#const mutation_stddev = 0.05
const mutation_stddev_list = [0.04]
#const num_attributes = 2        # number attributes for quantitative representation
const num_attributes_list = [2]        # number attributes for quantitative representation
const ngens = 3       # Generations after burn-in
const num_subpops = 1
const burn_in= 0    # if integer, generations of burn in.  If float, int_burn_in=burn_in*N+50
const ideal = 1.0
const fit_slope = 1.0
const neutral = false
