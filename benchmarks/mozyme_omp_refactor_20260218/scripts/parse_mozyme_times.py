#!/usr/bin/env python3
import argparse
import csv
import pathlib
import re
import statistics

PATTERNS = {
    "diagg1": re.compile(r"AFTER DIAGG1 IN ITE.*WALL-CLOCK INTERVAL:\s*([0-9.]+)"),
    "diagg2": re.compile(r"AFTER DIAGG2 IN ITE.*WALL-CLOCK INTERVAL:\s*([0-9.]+)"),
    "density": re.compile(r"After DENSIT\s+CPU INTERVAL:.*WALL-CLOCK INTERVAL:\s*([0-9.]+)"),
    "buildf": re.compile(r"After BUILDF\s+CPU INTERVAL:.*WALL-CLOCK INTERVAL:\s*([0-9.]+)"),
    "hof_kj": re.compile(r"FINAL HEAT OF FORMATION =.*=\s*(-?[0-9.]+) KJ/MOL"),
}


def parse_run(base: pathlib.Path):
    out_text = base.with_suffix(".out").read_text(errors="ignore")
    wall = float(base.with_suffix(".time").read_text().strip())

    vals = {}
    for key, rgx in PATTERNS.items():
        if key == "hof_kj":
            ms = list(rgx.finditer(out_text))
            vals[key] = float(ms[-1].group(1)) if ms else float("nan")
        else:
            vals[key] = sum(float(m.group(1)) for m in rgx.finditer(out_text))
    vals["wall"] = wall
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-dir", required=True)
    ap.add_argument("--threads", required=True, help="comma-separated, e.g. 1,2,4,8,16")
    ap.add_argument("--runs", type=int, default=3)
    args = ap.parse_args()

    bench = pathlib.Path(args.bench_dir)
    threads = [int(x) for x in args.threads.split(",") if x]

    rows = []
    for t in threads:
        for r in range(1, args.runs + 1):
            base = bench / f"receptor_diag2_autotune_t{t}_r{r}"
            vals = parse_run(base)
            rows.append({"threads": t, "run": r, **vals})

    raw = bench / "times_raw.csv"
    fields = ["threads", "run", "wall", "diagg1", "diagg2", "density", "buildf", "hof_kj"]
    with raw.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    metrics = ["wall", "diagg1", "diagg2", "density", "buildf", "hof_kj"]
    by = {}
    for row in rows:
        by.setdefault(row["threads"], []).append(row)

    summary = bench / "times_summary.csv"
    with summary.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["threads", "n"] + [f"{m}_mean" for m in metrics] + [f"{m}_std" for m in metrics])
        for t in sorted(by):
            g = by[t]
            vm = {m: [float(x[m]) for x in g] for m in metrics}
            w.writerow([t, len(g)] +
                       [f"{statistics.mean(vm[m]):.8f}" for m in metrics] +
                       [f"{statistics.stdev(vm[m]):.8f}" if len(vm[m]) > 1 else "0.00000000" for m in metrics])

    base = {m: statistics.mean([float(x[m]) for x in by[min(threads)]]) for m in metrics}
    scaling = bench / "times_scaling.csv"
    with scaling.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["threads", "wall_speedup", "diagg1_speedup", "diagg2_speedup", "density_speedup", "buildf_speedup"])
        for t in sorted(by):
            means = {m: statistics.mean([float(x[m]) for x in by[t]]) for m in metrics}
            w.writerow([
                t,
                f"{base['wall']/means['wall']:.8f}",
                f"{base['diagg1']/means['diagg1']:.8f}",
                f"{base['diagg2']/means['diagg2']:.8f}",
                f"{base['density']/means['density']:.8f}",
                f"{base['buildf']/means['buildf']:.8f}",
            ])

    print(raw)
    print(summary)
    print(scaling)


if __name__ == "__main__":
    main()
