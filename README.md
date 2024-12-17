# Matrix Inversion Optimization

This repository contains a C++ program to invert large dense matrices using Gauss-Jordan elimination. The project demonstrates performance optimizations including:
- Compiler optimizations (`-O3 -Ofast -march=native -ffast-math`)
- Semi-automatic vectorization (`#pragma clang loop vectorize(enable)`)
- Parallelization with OpenMP

**Steps to build and run:**
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./lab2
