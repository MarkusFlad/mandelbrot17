#!/bin/bash

# Compile the mandelbrot executable.
# First parameter: Name of the compiler.
# Second parameter: V1=no intrinsics version 1,
#                   V2=no intrinsics version 2,
#                   V3=no intrinsics version 3
#                   Default: use intrinsics
doCompile() {
    intrinsics=""
    executableName="$1-intrinsics"
    cppStd="c++20"
    if [[ $1 == "g++-9" ]]
      then
        cppStd="c++17"
    fi
    if [[ $2 == "V1" ]]
      then
        intrinsics=-DNO_INTRINSICS_V1
        executableName="$1-no-intrinsics-v1"
    elif [[ $2 == "V2" ]]
      then
        intrinsics=-DNO_INTRINSICS_V2
        executableName="$1-no-intrinsics-v2"
    fi
    set -x # echo on
    $1 --std=$cppStd -O3 -Wall -march=native -mno-fma $intrinsics src/mandelbrot17.cpp -lpthread -o Release/$executableName
    { set +x; } 2>/dev/null # echo off
}

# Run the mandelbrot performance test
# First parameter: Name of the compiler.
# Second parameter: V1=no intrinsics version 1,
#                   V2=no intrinsics version 2,
#                   V3=no intrinsics version 3
#                   Default: use intrinsics
runTest() {
    executableName="$1-intrinsics"
    if [[ $2 == "V1" ]]
      then
        executableName="$1-no-intrinsics-v1"
    elif [[ $2 == "V2" ]]
      then
        executableName="$1-no-intrinsics-v2"
    fi
    echo "Time for $executableName"
    python mandelbrot17.py $executableName
}

if [[ $1 != "NO_COMPILATION" ]]
    then
      doCompile g++-9
      doCompile g++-9 V1
      doCompile g++-9 V2
      doCompile g++-10
      doCompile g++-10 V1
      doCompile g++-10 V2
      doCompile g++-11
      doCompile g++-11 V1
      doCompile g++-11 V2
      doCompile clang++-10
      doCompile clang++-10 V1
      doCompile clang++-10 V2
      doCompile clang++-11
      doCompile clang++-11 V1
      doCompile clang++-11 V2
      doCompile clang++-12
      doCompile clang++-12 V1
      doCompile clang++-12 V2
fi

if [[ $1 != "NO_RUN" ]]
    then
      runTest g++-9
      runTest g++-9 V1
      runTest g++-9 V2
      runTest g++-10
      runTest g++-10 V1
      runTest g++-10 V2
      runTest g++-11
      runTest g++-11 V1
      runTest g++-11 V2
      runTest clang++-10
      runTest clang++-10 V1
      runTest clang++-10 V2
      runTest clang++-11
      runTest clang++-11 V1
      runTest clang++-11 V2
      runTest clang++-12
      runTest clang++-12 V1
      runTest clang++-12 V2
fi
