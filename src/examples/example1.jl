# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run_cv.jl examples/example1
or
>  julia -p 4 -L ContVarEvolution.jl run_cv.jl examples/example1
=#
export simtype
@everywhere simtype = 2    
const num_trials=1
#@everywhere const N = 8        # Meta-population size
@everywhere const N_list = [4]        # Meta-population size list
#const mutation_stddev = 0.05
#const num_attributes = 1        # number attributes for quantitative representation
#const num_attributes = 1        # number attributes for quantitative representation
const ngens = 2       # Generations after burn-in
const num_attributes_list = [2]        # number attributes for quantitative representation
const mutation_stddev_list = [0.05]
#N_mut_list = [1.0]
const num_subpops = 2
const burn_in= 0.1    # generations of burn_in as a multiple of N
const ideal = 1.0
#const mu = 0.0
const neutral = true
const wrap_attributes=false # wrap attribute values to stay within the interval [0,1]
const additive_error=false  # use additive error when mutating attributes as opposed to mulitiplicative error
