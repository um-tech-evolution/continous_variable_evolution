# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia hrun.jl examples/hexample1
or
>  julia -p 4 hrun.jl examples/hexample1
=#
@everywhere simtype = 4
const num_trials=3
#const N = 25          # Meta-population size
const N_list = [25,50,100,200]          # Meta-population size
const mutation_stddev = 0.080
const mutation_stddev_list = [0.080]
const N_mut_list = Float64[1.0]
const mutation_bias = 0.95
const num_attributes = 1        # number attributes for quantitative representation
const ngens = 20       # Generations after burn-in
const num_subpops = 1
const burn_in= 0    # if integer, generations of burn in.  If float, int_burn_in=burn_in*N+50
const ideal = 1.0
const fit_power = 3.0
const renormalize = true
const neutral = false
