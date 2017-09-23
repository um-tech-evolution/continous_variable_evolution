# Configuration for running spatial simulation
# simtype=2 means spatial structure by changing the ideal values for attributes
@everywhere simtype = 2    
@everywhere const N_list = [25,50,100]        # Meta-population size list
const num_trials = 1
const mutation_stddev_list = [0.005,0.01,0.02,0.04]
#const N_mut_list = [0.1,0.2,0.3 ]
const num_subpops=1
const num_attributes_list = [5]        # number attributes for quantitative representation
const ngens = 50       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const ideal = 0.5     # should be 0.5 for additive error
const neutral=false
const wrap_attributes=false # wrap attribute values to stay within the interval [0,1]
const additive_error=true
#const additive_error=false  # use additive error when mutating attributes as opposed to mulitiplicative error


