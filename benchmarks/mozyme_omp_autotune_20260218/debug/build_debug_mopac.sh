#!/usr/bin/env bash
set -euo pipefail

# Debug build flags requested in AGENTS.md
#   -g -Og -fopenmp -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
#   -finit-real=snan -finit-integer=-999999

build_dir=${1:-build-debug-diagg2}

cmake -S . -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_FLAGS_DEBUG='-g -Og -fopenmp -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999999'

cmake --build "$build_dir" -j"$(nproc)"

echo "Debug executable: $build_dir/mopac"
