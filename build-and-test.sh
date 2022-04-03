#!/bin/bash

# Compile the mandelbrot executable.
# First parameter: Name of the compiler.
doCompile() {
    preprocessorDefines=""
    executableName="mandelbrot-$1"
    set -x # echo on
    $1 -std=c++17 -O3 -Wall -march=native $preprocessorDefines src/mandelbrot17.cpp -lpthread -o Release/$executableName
    { set +x; } 2>/dev/null # echo off
}

# Run the mandelbrot performance test
# First parameter: Name of the compiler.
runTest() {
    executableName="mandelbrot-$1"
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
