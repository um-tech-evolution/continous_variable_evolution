# Configuration for running spatial simulation
export simtype
#global simtype=1 # means spatial structure by fitness or selection coefficient adjustment
# simtype=2 means spatial structure by changing the ideal values for attributes
@everywhere simtype = 2    
#@everywhere const N_list = [512,256,128,64,32,16,8,4,2]        # Meta-population size
#@everywhere const N_list = [128,64,32,16]        # Meta-population size
const num_trials = 1
@everywhere const N_list = [16]        # Meta-population size
#mutation_stddev_list = [0.005,0.01, 0.02, 0.05,0.1]
N_mut_list = [1.0, 2.0]
const num_subpops = 1                    # Number of subpopulations
const num_attributes_list = [1,5,10,50]        # number attributes for quantitative representation
const ngens = 200       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const ideal = 0.5
const wrap_attributes=false # wrap attribute values to stay within the interval [0,1]
const additive_error=false  # use additive error when mutating attributes as opposed to mulitiplicative error
const neutral = false



