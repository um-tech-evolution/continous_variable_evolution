# Configuration for running spatial simulation
export simtype
const num_trials = 2
simtype = 1    
const N_list = [4,8,16,32,64]        # Meta-population size list
const num_attributes_list = [5]        # number attributes for quantitative representation
const N_mut_list = [0.5, 1.0, 2.0, 4.0]
const num_subpops = 1
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const fit_slope = 1.0
const neutral = true


const w = 0.1
const a = 0.0
const b = 2.0
