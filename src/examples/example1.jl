# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run_cv.jl examples/example1
or
>  julia -p 4 -L ContVarEvolution.jl run_cv.jl examples/example1
=#
export simtype
@everywhere simtype = 2    
const num_trials=2
#@everywhere const N = 8        # Meta-population size
@everywhere const N_list = [10]        # Meta-population size list
#const mutation_stddev = 0.05
#const num_attributes = 1        # number attributes for quantitative representation
const ngens = 1       # Generations after burn-in
const num_attributes_list = [1]        # number attributes for quantitative representation
const mutation_stddev_list = [0.04]
#N_mut_list = [1.0]
const num_subpops = 1
const burn_in= 3    # if integer, generations of burn in.  If float, int_burn_in=burn_in*N+50
const ideal = 1.0
const fit_slope = 1.0
#const mu = 0.0
const neutral = true
const wrap_attributes=false # wrap attribute values to stay within the interval [0,1]
const additive_error=false  # use additive error when mutating attributes as opposed to mulitiplicative error
