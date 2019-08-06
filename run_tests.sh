#!/bin/bash
# Bash scrpt to run basic tests
cd src
julia run.jl examples/ent_example1 1
julia run.jl examples/ent_example2 1
julia run.jl examples/ent_example3 1
julia -p 4 run.jl examples/ent_example4 1
julia run.jl examples/ent_example5 1
julia run.jl examples/example1 1
julia run.jl examples/example2 1
julia run.jl examples/example3 1
julia -p 4 run.jl examples/example4 1
julia run.jl examples/example5 1
cd ..

