# MOZYME refactor benchmark summary

- New branch wall mean (t1): 56.868 s
- New branch wall mean (t16): 14.114 s
- Speedup t16/t1 (new): 4.029
- Speedup vs old baseline t1 at t16: 4.043

- Kernel speedup at t16 vs t1 (new):
  - diagg1: 4.609
  - diagg2: 2.580
  - density: 4.403
  - buildf: 5.631

- Kernel wall-fraction (mean):
  - t1: diagg1=0.394, diagg2=0.108, density=0.283, buildf=0.187
  - t16: diagg1=0.344, diagg2=0.169, density=0.259, buildf=0.134

- HOF spread across threads (mean values):
  - max-min = 0.032212 kJ/mol
