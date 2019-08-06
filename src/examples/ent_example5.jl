# Configuration for running spatial simulation
export simtype
#global simtype=1 # means spatial structure by fitness or selection coefficient adjustment
simtype = 1    
const num_trials = 1
const N_list = [16]        # Meta-population size
#mutation_stddev_list = [0.005,0.01, 0.02, 0.05,0.1]
N_mut_list = [1.0, 2.0]
const num_subpops = 1                    # Number of subpopulations
const num_attributes_list = [1,5,10,50]        # number attributes for quantitative representation
const ngens = 200       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const fit_slope = 1.0
const neutral = false



const w = 0.1
const a = 0.0
const b = 2.0
