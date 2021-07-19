#!/bin/bash

# Compile the mandelbrot executable.
# First parameter: Name of the compiler.
# Second parameter: No for using no intrinsics, else yes.
# Third parameter: Name of the executable
doCompile() {
    intrinsics=""
    if [[ $2 == "no" ]]
      then
        intrinsics=-DNO_INTRINSICS
    fi
    set -x # echo on
    $1 --std=c++17 -O3 -Wall -march=native -mno-fma $intrinsics src/mandelbrot17.cpp -lpthread -o Release/$3
    { set +x; } 2>/dev/null # echo off
}

# Run the mandelbrot performance test
# First parameter: The name of the executable
runTest() {
    echo "Time for $1"
    python mandelbrot17.py $1
}

doCompile g++-9 no gcc-9-no-intrinsics
doCompile g++-9 yes gcc-9-intrinsics
doCompile g++-10 no gcc-10-no-intrinsics
doCompile g++-10 yes gcc-10-intrinsics
doCompile g++-11 no gcc-11-no-intrinsics
doCompile g++-11 yes gcc-11-intrinsics

runTest gcc-9-no-intrinsics
runTest gcc-9-intrinsics
runTest gcc-10-no-intrinsics
runTest gcc-10-intrinsics
runTest gcc-11-no-intrinsics
runTest gcc-11-intrinsics
