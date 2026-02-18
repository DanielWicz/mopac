# perf commands for MOZYME/DIAGG2 analysis

Run these on the host shell (outside sandbox) when needed:

```bash
sudo sysctl -w kernel.perf_event_paranoid=0
sudo sysctl -w kernel.kptr_restrict=0
```

## Build with frame pointers
```bash
cmake -S . -B build-perf \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE='-O2 -g -fopenmp -fno-omit-frame-pointer'
cmake --build build-perf -j$(nproc)
```

## perf stat sweep
```bash
for t in 1 2 4 8 16; do
  OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=1 \
  perf stat -d -r 3 -- \
    ./build-perf/mopac receptor_mozyme_pm6_1scf_profile.mop
 done
```

## perf record + report
```bash
OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 \
perf record -g -- ./build-perf/mopac receptor_mozyme_pm6_1scf_profile.mop
perf report --stdio | head -n 200
```
