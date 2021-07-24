# mandelbrot17
Mandelbrot implementation in C++ that may use C++ features up to C++17.

This project produced the currently fastest version of the mandelbrot program called "C++g++#4" on the "Computer Language Benchmarks Game".
The current results for the various implementations of the mandelbrot program can be found under the following link:
https://benchmarksgame-team.pages.debian.net/benchmarksgame/performance/mandelbrot.html

To build with gcc you can just can follow the following steps:
```
    mkdir Release
    g++ --std=c++17 -O3 -Wall -march=native -mno-fma src/mandelbrot17.cpp -lpthread -o Release/mandelbrot17
```
    
And with the following command you can run the performance test that was basically done for the "Computer Language Benchmarks Game":
```
    python mandelbrot17.py
```
    
Recently I found that with a slight modification to the code, almost the same speed can be achieved without using SIMD intrinsics.
However, this only applies to Clang at the moment. So I created the script build-and-test.sh, which uses different versions of compilers
with and without intrinsics. The first version without intrinsics is the version that was uploaded to "Computer Language Benchmarks Game".
The second version is the one that is very fast in Clang. Note that gcc-11 also gives some better results with the second version, but is
still half as fast as Clang.

You can reproduce these tests by running the following lines on the command line (assuming folder Release exists):
```
    ./build-and-test.sh
```
    
This command creates and runs the tests. Once built, you can just use the following command to run the tests without building:
```
    ./build-and-test.sh NO_COMPILATION
```

On my machine I get the following results:
```
    Time for g++-9-intrinsics
    0.55745100975
    Time for g++-9-no-intrinsics-v1
    1.94521594048
    Time for g++-9-no-intrinsics-v2
    1.95812702179
    Time for g++-10-intrinsics
    0.556442975998
    Time for g++-10-no-intrinsics-v1
    1.95684885979
    Time for g++-10-no-intrinsics-v2
    1.96461892128
    Time for g++-11-intrinsics
    0.548613071442
    Time for g++-11-no-intrinsics-v1
    1.83499288559
    Time for g++-11-no-intrinsics-v2
    1.09055900574
    Time for clang++-10-intrinsics
    0.553403139114
    Time for clang++-10-no-intrinsics-v1
    2.5211930275
    Time for clang++-10-no-intrinsics-v2
    0.609521150589
    Time for clang++-11-intrinsics
    0.552711963654
    Time for clang++-11-no-intrinsics-v1
    2.5357940197
    Time for clang++-11-no-intrinsics-v2
    0.603178024292
    Time for clang++-12-intrinsics
    0.554489850998
    Time for clang++-12-no-intrinsics-v1
    2.53116297722
    Time for clang++-12-no-intrinsics-v2
    0.601706981659
```
