# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia run.jl examples/example1
or
>  julia -p 4 run.jl examples/example1
=#
export simtype
@everywhere simtype = 2    
const num_trials = 1
@everywhere const N_list = [50]        # Meta-population size
N_mut_list = [1.0]
#mutation_stddev_list = [0.005,0.01,0.02,0.04]
const num_subpops = 1                    # Number of subpopulations
const num_attributes_list = [1]        # number attributes for quantitative representation
const ngens = 1       # Generations after burn-in
const burn_in= 3.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const neutral = true
const fit_slope=1.0
const additive_error=false
const wrap_attributes=false
