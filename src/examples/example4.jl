# Configuration for running spatial simulation
export simtype
#global simtype=1 # means spatial structure by fitness or selection coefficient adjustment
# simtype=2 means spatial structure by changing the ideal values for attributes
const num_trials = 2
@everywhere simtype = 2    
@everywhere const N_list = [4,8,16,32,64]        # Meta-population size list
const num_attributes_list = [5]        # number attributes for quantitative representation
const N_mut_list = [0.5, 1.0, 2.0, 4.0]
const num_subpops = 1
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const fit_slope = 1.0
const neutral = true
const wrap_attributes=false # wrap attribute values to stay within the interval [0,1]
const additive_error=false  # use additive error when mutating attributes as opposed to mulitiplicative error


