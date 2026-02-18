# perf stat summary (receptor MOZYME PM6 1SCF)

- t1: elapsed=57.198s, IPC=2.822, cache_miss_rate=0.404, cs/s=8.8, mig/s=0.0
- t16: elapsed=15.242s, IPC=1.524, cache_miss_rate=0.351, cs/s=7775.4, mig/s=2052.5

Interpretation:
- IPC drop and high context-switch/migration rates at t16 indicate scheduling/runtime overhead and memory pressure are dominant beyond core utilization growth.
