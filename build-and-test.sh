#!/bin/bash

# Compile the mandelbrot executable.
# First parameter: Name of the compiler.
# Second parameter: V1=no intrinsics version 1,
#                   V2=no intrinsics version 2,
#                   V3=no intrinsics version 3
#                   Default: use intrinsics
doCompile() {
    preprocessorDefines=""
    executableName="mandelbrot-$1"
    set -x # echo on
    $1 --std=c++17 -O3 -Wall -march=native -mno-fma $preprocessorDefines src/mandelbrot17.cpp -lpthread -o Release/$executableName
    { set +x; } 2>/dev/null # echo off
}

# Run the mandelbrot performance test
# First parameter: Name of the compiler.
# Second parameter: V1=no intrinsics version 1,
#                   V2=no intrinsics version 2,
#                   V3=no intrinsics version 3
#                   Default: use intrinsics
runTest() {
    executableName="mandelbrot-$1"
    if [[ $2 == "StrangeSimdHint" ]]
      then
        executableName="mandelbrot-$1-strangeSimdHint"
    fi
    echo "Time for $executableName"
    python mandelbrot17.py $executableName
}

if [[ $1 != "NO_COMPILATION" ]]
    then
      doCompile g++-9
      doCompile g++-10
      doCompile g++-11
      doCompile clang++-10
      doCompile clang++-11
      doCompile clang++-12
fi

if [[ $1 != "NO_RUN" ]]
    then
      runTest g++-9
      runTest g++-10
      runTest g++-11
      runTest clang++-10
      runTest clang++-11
      runTest clang++-12
fi
