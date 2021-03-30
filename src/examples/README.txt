To run examples, do the following on the command line when in the src/ subdirectory:

> julia run.jl examples/example?

where ? is one of 1, 2, 3, 4, or 5.

To specify a random number seed, add it as an extra command-line parameter:

> julia run.jl examples/example? 8   # set random number seed to 8

Multiple processes can be used with the julia "-p" option:
(This is transparent on a multicore machine.)

> julia -p 4 run.jl examples/example?     # run with 4 processors (4 cores)

Copied from continuous_variable_evolution/README.md
