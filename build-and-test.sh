#!/bin/bash

# Compile the mandelbrot executable.
# First parameter: Name of the compiler.
# Second parameter: No for using no intrinsics, else yes.
# Third parameter: Name of the executable
doCompile() {
    intrinsics=""
    executableName="$1-intrinsics"
    if [[ $2 == "no" ]]
      then
        intrinsics=-DNO_INTRINSICS
        executableName="$1-no-intrinsics"
    fi
    set -x # echo on
    $1 --std=c++17 -O3 -Wall -march=native -mno-fma $intrinsics src/mandelbrot17.cpp -lpthread -o Release/$executableName
    { set +x; } 2>/dev/null # echo off
}

# Run the mandelbrot performance test
# First parameter: The name of the executable
runTest() {
    echo "Time for $1"
    python mandelbrot17.py $1
}

doCompile g++-9 no
doCompile g++-9 yes
doCompile g++-10 no
doCompile g++-10 yes
doCompile g++-11 no
doCompile g++-11 yes
doCompile clang++-10 no
doCompile clang++-10 yes
doCompile clang++-11 no
doCompile clang++-11 no
doCompile clang++-12 yes
doCompile clang++-12 yes

runTest g++-9-no-intrinsics
runTest g++-9-intrinsics
runTest g++-10-no-intrinsics
runTest g++-10-intrinsics
runTest g++-11-no-intrinsics
runTest g++-11-intrinsics
runTest g++-10-no-intrinsics
runTest g++-10-intrinsics
runTest g++-11-no-intrinsics
runTest g++-11-intrinsics
runTest clang++-10-no-intrinsics
runTest clang++-10-intrinsics
runTest clang++-11-no-intrinsics
runTest clang++-11-no-intrinsics
runTest clang++-12-intrinsics
runTest clang++-12-intrinsics
