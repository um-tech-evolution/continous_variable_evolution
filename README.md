Simulation of the evolution of continuous-valued (real-valued) traits for a paper on nearly neutral evolution.
There is one trait which can have multiple attributes.  Each attribute has a fixed ideal value, and fitness
is based on how close the attribute is to the ideal value.  Mutation corresponds to copy error as described in 
(Eerkins & Lipo 2005), (Hamilton & Buchanan 2009), (Rorabaugh 2014).  The code has provisions for multiple
subpopulations, but currently there is always only 1 subpopulation.

To run examples, do the following on the command line when in the src/ subdirectory:

> julia run.jl examples/example?

where ? is one of 1, 2, 3, 4, or 5.  

To run entropy examples, do the following on the command line when in the src/ subdirectory:

> julia run.jl examples/ent_example?

where ? is one of 1, 2, 3, 4, or 5.  

To specify a random number seed, add it as an extra command-line parameter:

> julia run.jl examples/example? 8   # set random number seed to 8

Multiple processes can be used with the julia "-p" option:
(This is transparent on a multicore machine.)

> julia -p 4 run.jl examples/example?     # run with 4 processors (4 cores)

Do "bash run_tests.sh" to run all examples in a Linux environment.
