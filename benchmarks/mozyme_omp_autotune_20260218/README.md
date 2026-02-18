# MOZYME OMP DIAGG2 Autotuning Benchmark Pack (2026-02-18)

This directory stores reproducible artifacts for the DIAGG2 adaptive scheduler work.

## Scope
- `src/MOZYME/diagg2.F90`: runtime autotuner for mode selection
  - level-sweep OpenMP loop vs dependency-task batching
- `src/MOZYME/diagg1.F90`: comments/documentation of dynamic work balancing and deterministic merge equations

## Autotuner model (DIAGG2)
Per thread-slot `T = diagg2_threads`, maintain online performance state:

1. **Normalized cost** per active pair
   - `c = elapsed / nactive`
2. **EWMA update**
   - `c_new = (1-α) * c_old + α * c_obs`
3. **Mode decision with hysteresis**
   - choose task mode when `c_task <= m * c_level`
4. **Batch adaptation**
   - `b_{k+1} = clamp(b_k + dir * step, b_min, b_max)`

Constants used in code:
- `α = 0.25`
- hysteresis `m = 1.02`
- exploration period = 24 calls (after both modes have at least 3 samples)

## Files
- `results/`:
  - `autotune_times_raw.csv`
  - `autotune_times_summary.csv`
  - `autotune_times_scaling.csv`
  - `autotune_vs_prev_task.csv`
  - `autotune_vs_old.csv`
  - `baseline_prev_task_summary.csv`
  - `baseline_old_summary.csv`
  - `autotune_t12_summary.txt`
- `scripts/`:
  - `run_mozyme_scaling.sh`
  - `parse_mozyme_times.py`
- `debug/`:
  - `build_debug_mopac.sh`
  - `perf_commands.md`

## Notes
- This pack contains benchmark summaries and reproducibility scripts only (not full raw `.out/.arc` run dumps).
- `AGENTS.md` is intentionally excluded from commit/push.
