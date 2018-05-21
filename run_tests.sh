#!/bin/bash
# Bash scrpt to run basic tests
cd src
julia run.jl examples/example1
julia run.jl examples/example2
julia run.jl examples/example3
julia -p 4 run.jl examples/example4
julia run.jl examples/example5
cd ..

