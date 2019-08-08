# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run_cv.jl examples/example1
or
>  julia -p 4 -L ContVarEvolution.jl run_cv.jl examples/example1
=#
export simtype
simtype = 1
#const N = 8        # Meta-population size
#const N_list = [25,50,100,200]        # Meta-population size list
const num_trials = 1
const N_list = [20]        # Meta-population size list
const num_attributes_list = [1]        # number attributes for quantitative representation
#const N_mut_list = [0.5, 1.0, 2.0, 4.0]
const N_mut_list = [2.0]
const num_subpops = 1
const ngens = 120      # Generations after burn-in
const burn_in= 0.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const fit_slope = 1.0
const neutral=false
const w = 0.1
