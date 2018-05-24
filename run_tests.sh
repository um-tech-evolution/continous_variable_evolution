#!/bin/bash
# Bash scrpt to run basic tests
cd src
/usr/bin/julia run.jl examples/example1 1
/usr/bin/julia run.jl examples/example2 1
/usr/bin/julia run.jl examples/example3 1
/usr/bin/julia -p 4 run.jl examples/example4 1
/usr/bin/julia run.jl examples/example5 1
cd ..

