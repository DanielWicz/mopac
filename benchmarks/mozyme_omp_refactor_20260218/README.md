# MOZYME OMP refactor benchmark pack (2026-02-18)

This pack captures benchmarks for the DIAGG1 sparse-candidate refactor and DIAGG2 schedule/scratch refactors.

## Reproduce scaling benchmark

```bash
bench_root=benchmarks/mozyme_omp_refactor_20260218
bdir=<fresh_build_dir>
$bench_root/scripts/run_mozyme_scaling.sh "$bdir/mopac" receptor_mozyme_pm6_1scf_profile.mop "$bench_root/results" 1,2,4,8,16 5
$bench_root/scripts/parse_mozyme_times.py --bench-dir "$bench_root/results" --threads 1,2,4,8,16 --runs 5
```

## Reproduce perf stat snapshots

Before `perf` runs, set:

```bash
sudo sysctl -w kernel.perf_event_paranoid=0
```

Then run:

```bash
bdir=<fresh_build_dir>
for t in 1 16; do
  base=benchmarks/mozyme_omp_refactor_20260218/debug/perf/receptor_perf_t${t}
  cp receptor_mozyme_pm6_1scf_profile.mop ${base}.mop
  sed -i -E "1s/THREADS=[0-9]+/THREADS=${t}/" ${base}.mop
  perf stat -e task-clock,cycles,instructions,cache-references,cache-misses,context-switches,cpu-migrations \
    -o ${base}.perfstat.txt \
    env OMP_NUM_THREADS=${t} OPENBLAS_NUM_THREADS=1 "$bdir/mopac" ${base}.mop > ${base}.stdout 2>&1
done
```

## Main outputs

- `results/times_summary.csv`
- `results/times_scaling.csv`
- `results/scaling_vs_old_t1.csv`
- `results/benchmark_summary.md`
- `debug/perf/perf_summary.md`
