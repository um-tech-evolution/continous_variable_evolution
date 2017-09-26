#!/bin/bash
# Bash scrpt to run basic tests
cd src
julia -L ContVarEvolution.jl run.jl examples/example1
julia -L ContVarEvolution.jl run.jl examples/example2
julia -L ContVarEvolution.jl run.jl examples/example3
julia -p 4 -L ContVarEvolution.jl run.jl examples/example4
julia -L ContVarEvolution.jl run.jl examples/example5
cd ..

