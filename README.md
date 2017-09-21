Simulation of the evolution of continuous-valued (real-valued) traits for a paper on nearly neutral evolution.
There is one trait which can have multiple attributes.  Each attribute has a fixed ideal value, and fitness
is based on how close the attribute is to the ideal value.  Mutation corresponds to copy error as described in 
(Eerkins & Lipo 2005), (Hamilton & Buchanan 2009), (Rorabaugh 2014).  The code has provisions for multiple
subpopulations, but currently there is always only 1 subpopulation.

As of 9/21/17, this repository requires julia version 6.0.

To run examples, do the following on the command line:

> julia -L ContVarEvolution.jl run.jl examples/example?

where ? is one of 1, 2, 3, 4, or 5.  

Multiple processes can be used with the julia "-p" option:

> julia -p 4 -L ContVarEvolution.jl run.jl examples/example?
