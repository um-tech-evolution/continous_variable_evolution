#!/bin/bash
# Bash scrpt to run basic tests
cd src
julia5 -L ContVarEvolution.jl run.jl examples/example1
julia5 -L ContVarEvolution.jl run.jl examples/example2
julia5 -L ContVarEvolution.jl run.jl examples/example3
julia5 -p 4 -L ContVarEvolution.jl run.jl examples/example4
julia5 -L ContVarEvolution.jl run.jl examples/example5
cd ..

